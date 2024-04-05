use std::thread::{self, ScopedJoinHandle};

use anyhow::Context;
use crossbeam_channel::{bounded, unbounded, Receiver, Sender, TryRecvError};
use rand::{distributions::Uniform, prelude::*};
use rand_pcg::Pcg64Mcg;

use rs_htslib::{
    hts::{HtsFile, HtsMode, HtsRead},
    sam::{
        parse_mod_tags, BamRec, CigarOp, MMParse, ModIter, SamHeader, SamReader, SeqBase, BAM_FDUP,
        BAM_FQCFAIL, BAM_FSECONDARY, BAM_FSUPPLEMENTARY, BAM_FUNMAP,
    },
};

use super::{brec_block::*, config::Config, output, reference::Reference};

/// The counts are:
///   0 - not modified
///   1 - 5mC or 5hmC
///   2 - 5mC
///   3 - 5hmC
pub(super) struct CpG {
    pub(super) offset: u32,
    pub(super) fwd_counts: [u32; 4],
    pub(super) rev_counts: [u32; 4],
}

impl CpG {
    fn new(offset: u32) -> Self {
        Self {
            offset,
            fwd_counts: [0; 4],
            rev_counts: [0; 4],
        }
    }

    fn add_counts(&mut self, other: &Self) {
        for (p1, p2) in self.fwd_counts.iter_mut().zip(self.rev_counts.iter()) {
            *p1 += *p2
        }
    }
}

pub(super) struct CountBlock {
    pub(super) idx: usize,
    pub(super) tid: usize,
    pub(super) start: usize,
    pub(super) cpg_sites: Vec<CpG>,
}

impl CountBlock {
    /// Initialize a new count block.  Make a list of CpG sites in the reference sequence.   Note
    /// that rf has the reference only for the region covered by the block.
    fn new(idx: usize, tid: usize, start: usize, rf: &[u8]) -> Self {
        let mut cpg_sites = Vec::new();

        // Unlikely to happen, but just in case...
        assert!(rf.len() < u32::MAX as usize, "Region too large");

        // Go through reference to find CpG sites
        for (i, c) in rf.windows(2).enumerate() {
            if c[0] == b'C' && c[1] == b'G' {
                cpg_sites.push(CpG::new(i as u32))
            }
        }
        Self {
            idx,
            tid,
            start,
            cpg_sites,
        }
    }

    pub(crate) fn find_site(&mut self, x: u32, reverse: bool) -> Option<&mut [u32]> {
        match self.cpg_sites.binary_search_by_key(&x, |c| c.offset) {
            Ok(i) => Some(if reverse {
                &mut self.cpg_sites[i].rev_counts
            } else {
                &mut self.cpg_sites[i].fwd_counts
            }),
            _ => None,
        }
    }

    pub(crate) fn end(&self) -> usize {
        self.start + self.cpg_sites.last().map_or(0, |c| c.offset as usize)
    }

    pub(crate) fn add_counts(&mut self, other: &mut Self) {
        assert!(self.start >= other.start);
        // Only add counts if overlapping
        if self.start <= other.end() && !self.cpg_sites.is_empty() {
            let delta = (self.start - other.start) as u32;
            let first_x = self.cpg_sites[0].offset + delta;
            let ix = other
                .cpg_sites
                .binary_search_by_key(&first_x, |c| c.offset)
                .expect("CpG site not found!");
            for (c1, c2) in other.cpg_sites[ix..].iter().zip(self.cpg_sites.iter_mut()) {
                c2.add_counts(c1)
            }
            if other.end() > self.end() {
                let last_x = self.cpg_sites.last().map_or(0, |c| c.offset) + delta;
                let ix = other
                    .cpg_sites
                    .binary_search_by_key(&last_x, |c| c.offset)
                    .expect("CpG site not found!");
                for mut c in other.cpg_sites.drain(ix + 1..) {
                    c.offset -= delta;
                    self.cpg_sites.push(c)
                }
            }
        }
    }
}

// Combine mod probabilities for the same position (i.e., for 5hmC and 5mC)
struct MethItr<'a, 'b> {
    mod_iter: ModIter<'a, 'b>,
    seq_ix: usize,
    finished: bool,
}

impl<'a, 'b> MethItr<'a, 'b> {
    fn new(mod_iter: ModIter<'a, 'b>) -> Self {
        Self {
            mod_iter,
            seq_ix: 0,
            finished: false,
        }
    }
}

impl<'a, 'b> Iterator for MethItr<'a, 'b> {
    type Item = (usize, SeqBase, Option<(u8, u8, bool)>);

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            None
        } else if let Some(d) = self.mod_iter.next_pos() {
            let i = self.seq_ix;
            let b = d.seq_base();
            let v = d.data();
            let p = if v.is_empty() {
                None
            } else {
                let (tp_m, tp_h) =
                    v.iter()
                        .fold((0, 0), |(t_m, t_h), (p, m)| match m.base_mod_code() {
                            Some(b'm') => (t_m + (*p as usize), t_h),
                            Some(b'h') => (t_m, t_h + (*p as usize)),
                            _ => (t_m, t_h),
                        });
                Some((
                    tp_m.min(255) as u8,
                    tp_h.min(255) as u8,
                    v[0].1.reverse_strand(),
                ))
            };
            self.seq_ix += 1;
            Some((i, b, p))
        } else {
            self.finished = true;
            None
        }
    }
}

fn get_read_coords(rec: &BamRec) -> anyhow::Result<(usize, usize, usize)> {
    let tid = rec.tid().ok_or(anyhow!("Missing tid for BAM record"))?;
    // If the previous command did not fail, then these should also work
    let x = rec.pos().unwrap() as usize;
    let y = rec.endpos().unwrap() as usize;
    Ok((tid, x, y))
}

fn process_record(
    cfg: &Config,
    rf: &[u8],
    rec: &BamRec,
    mm_parse: &mut MMParse,
    cb: &mut CountBlock,
) -> anyhow::Result<()> {
    // Find and parse MM/ML tags from record
    if let Some(tags) = parse_mod_tags(rec, mm_parse)? {
        // Tags found.  Make iterator over modified positions
        let mod_iter = mm_parse.mk_pos_iter(rec, &tags)?;
        let mut meth_itr = MethItr::new(mod_iter);

        // Threshold (on prob. scale of 0-255) for calling methylated or non-methylation sites
        let thresh = cfg.prob_threshold();
        let thresh1 = 255 - cfg.prob_threshold();
        // Sanity check
        assert!(thresh1 < thresh);

        // We go through the read following the CIGAR string, so we can match up the reference and
        // the read and match these up to the reported modifications

        // Get start of current count block
        let start_x = cb.start;

        // Get map coordinates from record (we should only have mapped reads here)
        let (_, x, y) = get_read_coords(rec)?;

        // Is read reversed in record compared to its original (sequenced) orientation
        let read_rev = rec.is_reversed();

        // Sanity check - we should only have reads that come at of after the start of the count block
        assert!(start_x <= x);

        // Difference between the starting map position and the start of the count block
        let delta_x = x - start_x;

        // Get reference seq and make an iterator over it
        let r = &rf[delta_x..y - start_x];
        let mut ref_itr = r.iter().enumerate();

        // Iterate through cigar ops
        let cigar = rec.cigar().ok_or(anyhow!("No CIGAR for read"))?;
        for c in cigar.iter() {
            let (op, l) = c.op_pair();

            // Not sure why we would have zero length cigar ops, but if we do we will ignore them
            if l == 0 {
                continue;
            }

            match op {
                // Matching ops
                CigarOp::Match | CigarOp::Diff | CigarOp::Equal => {
                    for _ in 0..l {
                        // Get next base from read along with any modification if present.
                        // If None is returned from meth_itr.next() then there are no more
                        // modifications present for this record.
                        if let Some((_, _, p)) = meth_itr.next() {
                            let (i, _) = ref_itr.next().expect("Bad CIGAR");
                            if let Some((m, rev)) = p.and_then(|(x_m, x_h, rev)| {
                                if x_m >= thresh {
                                    Some((2, rev))
                                } else if x_h >= thresh {
                                    Some((3, rev))
                                } else if x_m + x_h >= thresh {
                                    Some((1, rev))
                                } else if x_m + x_h <= thresh1 {
                                    Some((0, rev))
                                } else {
                                    None
                                }
                            }) {
                                // XOR
                                let r = (read_rev || rev) && !(read_rev && rev);

                                if let Some(z) = if r {
                                    if i > 0 {
                                        Some(i + delta_x - 1)
                                    } else {
                                        None
                                    }
                                } else {
                                    Some(i + delta_x)
                                } {
                                    if let Some(ct) = cb.find_site(z as u32, r) {
                                        ct[m] += 1;
                                        if m >= 2 {
                                            ct[1] += 1;
                                        }
                                    }
                                }
                            }
                        } else {
                            // No more modified positions for this record
                            break;
                        }
                    }
                }
                _ => {
                    let t = c.op_type();
                    if (t & 1) == 1 {
                        // Consumes query
                        if meth_itr.nth((l - 1) as usize).is_none() {
                            break;
                        }
                    }
                    if (t & 2) == 2 {
                        // Consumes reference
                        ref_itr.nth((l - 1) as usize);
                    }
                }
            }
        }
    }

    Ok(())
}

fn read_thread(
    cfg: &Config,
    ix: usize,
    contigs: &[String],
    reference: &Reference,
    out_send: Sender<CountBlock>,
    block_recv: Receiver<ProcessBlock>,
    block_send: Sender<BRecBlock>,
) -> anyhow::Result<()> {
    debug!("read thread {} starting up", ix);

    let mut mm_parse = MMParse::default();
    mm_parse
        .set_selection(&["C+m", "C+h"])
        .expect("Problem setting modification selection");

    for proc_blk in block_recv.iter() {
        let ProcessBlock {
            idx,
            end_pos,
            mut bblock,
        } = proc_blk;
        trace!("read thread {} received block {}", ix, idx);

        let count_block = process_brec_block(
            cfg,
            contigs,
            reference,
            &mut mm_parse,
            idx,
            &mut bblock,
            end_pos,
        )?;

        // Send count block to output thread
        out_send.send(count_block)?;
        // Send completed block to back to main thread
        block_send.send(bblock)?
    }
    debug!("read thread {} shutting down", ix);
    Ok(())
}

fn process_brec_block(
    cfg: &Config,
    contigs: &[String],
    reference: &Reference,
    mm_parse: &mut MMParse,
    idx: usize,
    bblock: &mut BRecBlock,
    end_pos: usize,
) -> anyhow::Result<CountBlock> {
    assert!(!bblock.is_empty());
    let recs = bblock.brec_vec();
    let (tid, start, _) = get_read_coords(&recs[0])?;
    assert!(start <= end_pos);
    let rf = reference
        .get_seq(contigs[tid].as_str(), start, end_pos)
        .ok_or(anyhow!(
            "Could not get reference sequence for contig {}",
            contigs[tid]
        ))?;

    trace!(
        "Setting up count block for {}:{}-{}",
        contigs[tid].as_str(),
        start,
        end_pos
    );
    let mut cb = CountBlock::new(idx, tid, start, rf);

    for rec in recs {
        process_record(cfg, rf, rec, mm_parse, &mut cb)?
    }

    Ok(cb)
}

/// The main event.  Read input (SAM/BAM/CRAM from file or stdin)
/// If a record is mapped and has MM and ML tags, extract methylation information and output
/// a bedmethyl (CpG) file.
///
/// Threading model
///
/// We have a pool of hts threads set by cfg.threads(), the main thread, an output thread
/// and a number of process threads.
///
/// The main thread will read BamRec from input and store in blocks of size BREC_BLOCK_SIZE, which
/// will be sent to the process threads to extract the methylation information and get counts
/// of methylated/non-methylated per CpG.  All reads in the block will map to the same contig.
/// Counts from completed blocks will be sent to the output thread which will assemble them and
/// output the Bedmethyl file. Finished blocks of BamRecs will be sent back to the main thread so
/// that they can be reused.
pub fn read_input(cfg: &Config, rf: &Reference) -> anyhow::Result<()> {
    debug!("Processing input");
    let mut hts_file =
        HtsFile::open(cfg.input_file(), HtsMode::Read).with_context(|| "Error opening input")?;
    trace!("Opened input successfully");
    if let Some(tpool) = cfg.hts_thread_pool() {
        hts_file.attach_thread_pool(tpool)?;
    }
    let hdr = SamHeader::read(&mut hts_file).with_context(|| "Error reading input header")?;
    trace!("Read in input header");

    // Collect contig names
    let seq = hdr.sequences();
    let mut ctg_names: Vec<_> = Vec::with_capacity(seq.size_hint().0);

    for (s, _) in seq {
        ctg_names.push(s.to_str()?.to_owned())
    }

    // Set number of process threads, so we have 1-2 hts threads per proc thread
    let n_proc = (cfg.threads() + 1) / 2;

    // Create list of blocks to hold records
    // We create enough so that each process thread and the output thread can have 1 block being
    // worked on and one in reserve
    let nb = (n_proc + 1) << 1;
    let mut b_rec_blocks = Vec::with_capacity(nb);
    for _ in 0..nb {
        b_rec_blocks.push(BRecBlock::new())
    }

    // Create channels

    // Channels for communication between the main thread and the process threads
    // These are for sending blocks to the process threads
    let (block_send, block_recv) = unbounded();
    // and these are for sending used blocks back to the main thread
    let (used_send, used_recv) = unbounded();
    // These channels are for sending count data from the process threads to the output thread
    let (out_send, out_recv) = bounded(n_proc * 8);

    // Create reader
    let mut rdr = SamReader::new(hts_file, hdr);

    // Create RNG for down sampling
    let mut rng = if cfg.discard() > 0.0 {
        Some(Pcg64Mcg::seed_from_u64(cfg.seed()))
    } else {
        None
    };

    let mut res = false;
    // Create threads
    thread::scope(|s| {
        // Spawn output thread
        let ctg_ref = &ctg_names;
        let output = s.spawn(|| output::output_thread(cfg, ctg_ref, out_recv));

        // Spawn process threads
        let mut process: Vec<_> = (0..n_proc)
            .map(|ix| {
                let sc = out_send.clone();
                let rc = block_recv.clone();
                let uc = used_send.clone();
                let ctg_ref = &ctg_names;

                s.spawn(move || read_thread(cfg, ix, ctg_ref, rf, sc, rc, uc))
            })
            .collect();

        // We do this so that when the process threads exit the channel will be disconnected so the
        // output thread will exit
        drop(out_send);

        // Consecutive index for ProcessBlock structs
        let mut output_idx = 0;

        // Main process loop.  Returns true if an error occurs
        res = loop {
            // Check if any of the process threads are finished.  If they have at this stage this
            // must be because they hit an error condition, so we abort at this point
            if process.iter_mut().any(|jh| jh.is_finished()) {
                break true;
            }

            // Get an empty BRecBlock to fill.  If None is returned than an error occurred, so we abort
            let Some(mut blk) = get_b_rec_block(&mut b_rec_blocks, &used_recv) else {
                break true;
            };

            // Fill BRecBlock with BamRec records from input stream
            let (eof, end_pos) = match fill_b_rec_block(cfg, &mut rdr, rng.as_mut(), &mut blk) {
                Ok(x) => x,
                Err(e) => {
                    error!("Error returned from reading BAM records: {}", e);
                    break true;
                }
            };

            if blk.is_empty() {
                assert!(eof, "Empty block");
                break false;
            }

            // Create ProcessBlock and send to process threads
            let pb = ProcessBlock::new(blk, output_idx, end_pos.expect("Missing end position"));
            trace!("Sending block {} to process threads", output_idx);
            output_idx += 1;
            if let Err(e) = block_send.send(pb) {
                error!("Error sending process block to process threads: {}", e);
                break true;
            }

            // If end of file/stream reached we terminate normally
            if eof {
                break false;
            }
        };

        // Drop block_send to signal to process threads that input is finished
        drop(block_send);

        // Wait until process threads have finished and recover any errors
        res = join_process_and_output_threads(process, output) || res;
    });

    if res {
        Err(anyhow!("Error - reading of input file/stream unsuccessful"))
    } else {
        info!("Finished processing input");
        Ok(())
    }
}

/// Fill BRecBlock blk with BamRec read from hts. Returns an indicator of end-of-stream
/// or an error if any error occurred while reading
fn fill_b_rec_block(
    cfg: &Config,
    rdr: &mut SamReader,
    mut rng: Option<&mut Pcg64Mcg>,
    blk: &mut BRecBlock,
) -> anyhow::Result<(bool, Option<usize>)> {
    // Reset block to empty state
    blk.clear();

    let mut curr: Option<(usize, usize)> = None;

    let res = loop {
        // Get next available BamRec or terminate normally
        let Some(rec) = blk.next_rec() else {
            break false;
        };

        // Read BamRec.
        if !rdr.read(rec)? {
            // End of file/stream
            // This is needed so that the failed read is not output
            blk.decr_ix();
            break true;
        }

        // If we are down sampling then we randomly discard reads *before* looking at them
        let discard = rng
            .as_mut()
            .map(|z| z.sample(Uniform::from(0.0..1.0)) < cfg.discard())
            .unwrap_or(false);

        // For reads that are not down sampled, we filter on the flags and on MAPQ
        if discard
            || rec.flag_check_any(
                BAM_FUNMAP | BAM_FDUP | BAM_FQCFAIL | BAM_FSECONDARY | BAM_FSUPPLEMENTARY,
            )
            || rec.qual().unwrap() < cfg.mapq_threshold()
        {
            // remove discarded read from block
            blk.decr_ix();
        } else {
            // Check we are still on the same contig
            let tid = rec.tid().expect("No tid for mapped record");
            let y = rec.endpos().unwrap() as usize;
            if let Some((i, y1)) = curr.as_ref() {
                if *i != tid {
                    blk.decr_ix();
                    break false;
                }
                curr = Some((tid, y.max(*y1)));
            } else {
                curr = Some((tid, y))
            }
        }
    };
    let end_pos = curr.map(|c| c.1);
    Ok((res, end_pos))
}

/// Returns an empty BRecBlock from b_rec_blocks and recovers empty blocks from the output thread.
/// Will block if no BRecBlock available, and will return None on error.
fn get_b_rec_block(
    b_rec_blocks: &mut Vec<BRecBlock>,
    r: &Receiver<BRecBlock>,
) -> Option<BRecBlock> {
    match _get_block(b_rec_blocks, r) {
        Ok(b) => Some(b),
        Err(e) => {
            error!("Error received when recovering used blocks: {}", e);
            None
        }
    }
}

/// Returns an empty BRecBlock from b_rec_blocks and recovers empty blocks from the output thread.
/// Will block if no BRecBlock is available, and will return an error if the channel is disconnected.
fn _get_block(
    b_rec_blocks: &mut Vec<BRecBlock>,
    r: &Receiver<BRecBlock>,
) -> anyhow::Result<BRecBlock> {
    // Non-blocking recovery of used BRecBlocks from r
    try_recover_used_blocks(b_rec_blocks, r)?;

    loop {
        // Return BRecBlock if available
        if let Some(b) = b_rec_blocks.pop() {
            return Ok(b);
        }

        // Blocking recovery of used BRecBLocks from r
        recover_used_block(b_rec_blocks, r)?;
    }
}

/// If any used BRecBlock have been sent from the output thread, recover from r and add to
/// brec_blocks.  Will recover all BRecBlocks that are available in the channel, but will not block.
/// Returns an error if the channel is disconnected.
fn try_recover_used_blocks(
    b_rec_blocks: &mut Vec<BRecBlock>,
    r: &Receiver<BRecBlock>,
) -> anyhow::Result<()> {
    loop {
        match r.try_recv() {
            Ok(b) => b_rec_blocks.push(b),
            Err(TryRecvError::Empty) => break,
            Err(_) => return Err(anyhow!("Error - used block channel disconnected")),
        }
    }
    Ok(())
}

/// Recover used BRecBlocks from r and add to brec_blocks. Will block until a BRecBlock is received.
/// Returns an error if the channel is disconnected or some other error occurs on recv()
fn recover_used_block(
    b_rec_blocks: &mut Vec<BRecBlock>,
    r: &Receiver<BRecBlock>,
) -> anyhow::Result<()> {
    let b = r.recv().with_context(|| "Error receiving used blocks")?;
    b_rec_blocks.push(b);
    Ok(())
}

fn join_thread<F>(h: ScopedJoinHandle<anyhow::Result<()>>, f: F) -> bool
where
    F: FnOnce() -> String,
{
    if let Ok(r) = h.join() {
        if let Err(e) = r {
            error!("{} returned an error: {}", f(), e);
            true
        } else {
            false
        }
    } else {
        error!("{} panicked", f());
        true
    }
}

/// Joins all process threads and recover errors.  Returns true if any thread returned an error.
fn join_process_and_output_threads(
    mut v: Vec<ScopedJoinHandle<anyhow::Result<()>>>,
    out: ScopedJoinHandle<anyhow::Result<()>>,
) -> bool {
    let mut res = false;

    // Note we must join the process threads before we join the output thread otherwise
    // we will have a deadlock
    for (ix, jh) in v.drain(..).enumerate() {
        res = join_thread(jh, || format!("process thread {}", ix)) || res
    }
    res = join_thread(out, || "Output thread".to_owned()) || res;
    res
}
