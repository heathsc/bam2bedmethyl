use std::thread::{self, ScopedJoinHandle};

use anyhow::Context;
use crossbeam_channel::{unbounded, Receiver, Sender, TryRecvError};
use rand::{distributions::Uniform, prelude::*};
use rand_pcg::Pcg64Mcg;

use rs_htslib::{
    hts::{HtsFile, HtsMode, HtsRead},
    sam::{
        find_mod_tags, parse_mod_tags, BamAux, BamAuxItem, BamRec, CigarOp, MMParse, ModIter,
        Modification, SamHeader, SamReader, SeqBase, BAM_FDUP, BAM_FQCFAIL, BAM_FSECONDARY,
        BAM_FSUPPLEMENTARY, BAM_FUNMAP,
    },
    HStr,
};

use super::{brec_block::*, config::Config, output, reference::Reference};

pub(super) struct CpG {
    offset: u32,
    counts: [u32; 2],
}

impl CpG {
    fn new(offset: u32) -> Self {
        Self {
            offset,
            counts: [0; 2],
        }
    }
}

pub(super) struct CountBlock {
    tid: usize,
    start: usize,
    cpg_sites: Vec<CpG>,
}

impl CountBlock {
    /// Initialize a new count block.  Make a list of CpG sites in the reference sequence.   Note
    /// that rf has the reference only for the region covered by the block.
    fn new(tid: usize, start: usize, rf: &[u8]) -> Self {
        let mut cpg_sites = Vec::new();

        // Unlikely to happen, but just in case...
        assert!(rf.len() < u32::MAX as usize, "Region to large");

        // Go through reference to find CpG sites
        for (i, c) in rf.windows(2).enumerate() {
            if c[0] == b'C' && c[1] == b'G' {
                cpg_sites.push(CpG::new(i as u32))
            }
        }
        Self {
            tid,
            start,
            cpg_sites,
        }
    }
}

// Combine mod probs for the same position (i.e., for 5hmC and 5mC)
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
    type Item = (usize, SeqBase, Option<u8>);

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            None
        } else {
            if let Some(d) = self.mod_iter.next_pos() {
                let i = self.seq_ix;
                let b = d.seq_base();
                let v = d.data();
                let p = if v.is_empty() {
                    None
                } else {
                    let tp = v.iter().fold(0, |t, (p, m)| match m.base_mod_code() {
                        Some(b'm') | Some(b'h') => t + (*p as usize),
                        _ => t,
                    });
                    Some(tp.min(255) as u8)
                };
                self.seq_ix += 1;
                Some((i, b, p))
            } else {
                self.finished = true;
                None
            }
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
    rf: &Reference,
    rec: &BamRec,
    mm_parse: &mut MMParse,
    ctg_names: &[String],
) -> anyhow::Result<()> {
    if let Some(tags) = parse_mod_tags(rec, mm_parse)? {
        let mod_iter = mm_parse.mk_pos_iter(rec, &tags)?;
        let mut meth_itr = MethItr::new(mod_iter);

        let (tid, x, y) = get_read_coords(rec)?;
        let ctg = ctg_names[tid].as_str();
        let r = rf.get_seq(ctg, x, y).ok_or(anyhow!(
            "Could not get sequence for {}:{}-{}",
            ctg,
            x,
            y
        ))?;

        let mut ref_itr = r.iter().enumerate();
        let cigar = rec.cigar().ok_or(anyhow!("No CIGAR for read"))?;

        println!("Read {}\t{}", rec.qname(), rec.is_reversed());

        for c in cigar.iter() {
            let (op, l) = c.op_pair();
            if l == 0 {
                continue;
            }
            match op {
                CigarOp::Match | CigarOp::Diff | CigarOp::Equal => {
                    for _ in 0..l {
                        if let Some((_, b, p)) = meth_itr.next() {
                            let (i, c) = ref_itr.next().expect("Bad CIGAR");
                            if p.is_some() {
                                println!("{}\t{}\t{}\t{}\t{:?}", ctg, x + i, *c as char, b, p);
                            }
                        } else {
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

    for mut proc_blk in block_recv.iter() {
        let ProcessBlock {
            idx,
            end_pos,
            mut bblock,
        } = proc_blk;
        trace!("read thread {} received block {}", ix, idx);

        let count_block =
            process_brec_block(cfg, contigs, reference, &mut mm_parse, &mut bblock, end_pos)?;

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
    bblock: &mut BRecBlock,
    end_pos: usize,
) -> anyhow::Result<CountBlock> {
    assert!(!bblock.is_empty());
    let recs = bblock.brec_vec();
    let (tid, start, _) = get_read_coords(&recs[0])?;
    assert!(start >= end_pos);
    let rf = reference
        .get_seq(contigs[tid].as_str(), start, end_pos)
        .ok_or(anyhow!("Could not get reference for read block"))?;

    let mut cb = CountBlock::new(tid, start, rf);

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
pub fn read_input(cfg: &Config, rf: &Reference, mm_parse: &mut MMParse) -> anyhow::Result<()> {
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

    // Set number of process threads so we have 1-2 hts threads per proc thread
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
    let (out_send, out_recv) = unbounded();

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

            // Get an empty BRecBlock to fill.  If None is returned than an error occurred so we abort
            let Some(mut blk) = get_b_rec_block(&mut b_rec_blocks, &used_recv) else {
                break true
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
            break false
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
            .unwrap_or(true);

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
            if let Some((i, y1)) = curr.take() {
                if i != tid {
                    blk.decr_ix();
                    break false;
                }
                curr = Some((tid, y1.max(y)))
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
    // Non blocking recovery of used BRecBlocks from r
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
