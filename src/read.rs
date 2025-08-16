use std::thread::{self, ScopedJoinHandle};

use anyhow::Context;
use crossbeam_channel::{bounded, unbounded, Receiver, TryRecvError};
use rand::{distributions::Uniform, prelude::*};
use rand_pcg::Pcg64Mcg;

use brec_block::*;
use rs_htslib::{
    hts::{HtsFile, HtsMode, HtsPos},
    sam::{
        SamHeader, SamReader, BAM_FDUP, BAM_FQCFAIL, BAM_FSECONDARY, BAM_FSUPPLEMENTARY, BAM_FUNMAP,
    },
};

use super::{
    config::Config,
    output,
    process_read::{process_block::ProcessBlock, process_read_thread},
    read::read_record::ReadRec,
    reference::Reference,
};

pub mod brec_block;
mod opt_index;
pub mod pileup;
pub mod read_record;

use pileup::Pileup;

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
    let n_proc = cfg.threads().div_ceil(2);

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

                s.spawn(move || process_read_thread(cfg, ix, ctg_ref, rf, sc, rc, uc))
            })
            .collect();

        // We do this so that when the process threads exit the channel will be disconnected so the
        // output thread will exit
        drop(out_send);

        // Consecutive index for ProcessBlock structs
        let mut output_idx = 0;

        // Handle pileup if required
        let mut pileup = if cfg.pileup() {
            Some(Pileup::new())
        } else {
            None
        };
        let mut pending = None;

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
            let (eof, end_pos) = match fill_b_rec_block(
                cfg,
                &mut rdr,
                rng.as_mut(),
                &mut blk,
                pileup.as_mut(),
                &mut pending,
                &ctg_names,
                rf,
            ) {
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
#[allow(clippy::too_many_arguments)]
fn fill_b_rec_block(
    cfg: &Config,
    rdr: &mut SamReader,
    mut rng: Option<&mut Pcg64Mcg>,
    blk: &mut BRecBlock,
    mut pileup: Option<&mut Pileup>,
    pending: &mut Option<ReadRec>,
    ctg_names: &[String],
    rf: &Reference,
) -> anyhow::Result<(bool, Option<usize>)> {
    // Reset block to empty state
    blk.clear();

    let mut curr: Option<(usize, HtsPos, &[u64])> = None;

    let res = loop {
        // Get next available BamRec or terminate normally
        let rec = match blk.next_rec(rdr, pending)? {
            BrecFill::EndOfBlock => break false,
            BrecFill::EndOfFile => break true,
            BrecFill::Rec(r) => r,
        };

        // If we are down sampling then we randomly discard reads *before* looking at them
        let discard = rng
            .as_mut()
            .map(|z| z.sample(Uniform::from(0.0..1.0)) < cfg.discard())
            .unwrap_or(false);

        let br = rec.brec();

        // For reads that are not down sampled, we filter on the flags and on MAPQ
        if discard
            || br.flag_check_any(
                BAM_FUNMAP | BAM_FDUP | BAM_FQCFAIL | BAM_FSECONDARY | BAM_FSUPPLEMENTARY,
            )
            || br.qual().unwrap() < cfg.mapq_threshold()
        {
            // remove discarded read from block
            blk.decr_ix();
        } else {
            // Check we are still on the same contig
            let tid = br.tid().expect("No tid for mapped record");
            let y = br.endpos().unwrap();

            let (old_tid, end_pos, cpgs) = curr.take().unwrap_or_else(|| {
                let cpgs = rf
                    .get_cpgs(ctg_names[tid].as_str())
                    .expect("Missing CpG list for contig");
                (tid, 0, cpgs)
            });

            if old_tid != tid {
                *pending = Some(blk.push_back());
                if let Some(p) = pileup {
                    p.clear()
                }
                curr = Some((tid, end_pos, cpgs));
                break false;
            }

            let pileup_ix = if let Some(pl) = pileup.as_mut() {
                pl.trim(br.pos().unwrap() as usize);
                let y1 = y as u64;
                let cpg_ix = match cpgs.binary_search(&y1) {
                    Ok(i) => i + 1,
                    Err(i) => i,
                }
                .min(cpgs.len() - 1);
                let y1 = (cpgs[cpg_ix] as HtsPos + 1).max(y);
                Some(pl.new_index(y1 as usize))
            } else {
                None
            };
            rec.update(y, pileup_ix);
            curr = Some((tid, (y + 1).max(end_pos), cpgs))
        }
    };
    let end_pos = curr.map(|(_, y, _)| y as usize);
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
