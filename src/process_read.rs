pub mod count_block;
pub mod cytosine;
mod meth_itr;
pub mod process_block;

use crossbeam_channel::{Receiver, Sender};

use rs_htslib::sam::{BamRec, CigarOp, MMParse, parse_mod_tags};

use super::{
    config::Config,
    read::{
        pileup::{PileupCode, PileupEntry},
        read_record::ReadRec,
    },
    reference::Reference,
};

use crate::read::brec_block::*;

use count_block::CountBlock;
use meth_itr::MethItr;
use process_block::ProcessBlock;

fn get_read_coords(rec: &BamRec) -> anyhow::Result<(usize, usize)> {
    let tid = rec.tid().ok_or(anyhow!("Missing tid for BAM record"))?;
    // If the previous command did not fail, then these should also work
    let x = rec.pos().unwrap() as usize;
    Ok((tid, x))
}

const PILEUP_ITEMS: [PileupCode; 8] = [
    PileupCode::CytosineFwd,
    PileupCode::TotalMethFwd,
    PileupCode::MethCytosineFwd,
    PileupCode::HydroxyMethCytosineFwd,
    PileupCode::CytosineRev,
    PileupCode::TotalMethRev,
    PileupCode::MethCytosineRev,
    PileupCode::HydroxyMethCytosineRev,
];

fn process_record(
    cfg: &Config,
    rf: &[u8],
    rrec: &ReadRec,
    mm_parse: &mut MMParse,
    cb: &mut CountBlock,
) -> anyhow::Result<()> {
    let pileup_ix = rrec.pileup_ix();
    let y = rrec.end_pos() as usize;
    let rec = rrec.brec();
    // Find and parse MM/ML tags from record
    if let Some(tags) = parse_mod_tags(rec, mm_parse)? {
        let mut cts: [u32; 4] = [0; 4];

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
        let start_x = cb.start();

        // Get map coordinates from record (we should only have mapped reads here)
        let x = rec.pos().unwrap() as usize;

        // Is read reversed in record compared to its original (sequenced) orientation
        let read_rev = rec.is_reversed();

        // Sanity check - we should only have reads that come at of after the start of the count block
        assert!(start_x <= x);

        // Difference between the starting map position and the start of the count block
        let delta_x = x - start_x;

        // Get reference seq and make an iterator over it
        let r = &rf[delta_x..y - start_x];
        let mut ref_itr = r.iter().enumerate();

        // Setup pileup counts if required
        if let Some(ix) = pileup_ix {
            cb.init_pileup_index(ix, x, y)
        }

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
                            if let Some((om, x_m, x_h, rev)) = p.map(|(x_m, x_h, rev)| {
                                if x_m >= thresh {
                                    (Some(2), x_m, x_h, rev)
                                } else if x_h >= thresh {
                                    (Some(3), x_m, x_h, rev)
                                } else if x_m + x_h >= thresh {
                                    (Some(1), x_m, x_h, rev)
                                } else if x_m + x_h <= thresh1 {
                                    (Some(0), x_m, x_h, rev)
                                } else {
                                    (None, x_m, x_h, rev)
                                }
                            }) {
                                // XOR
                                let r = (read_rev || rev) && !(read_rev && rev);

                                if let Some(z) = if r {
                                    if i + delta_x > 0 {
                                        Some(i + delta_x - 1)
                                    } else {
                                        None
                                    }
                                } else {
                                    Some(i + delta_x)
                                } && let Some(cpg) = cb.find_site(z as u32, r)
                                {
                                    if let Some(m) = om {
                                        cts[m] += 1;
                                        cpg.incr_count(m);
                                        if let Some(k) = pileup_ix {
                                            cpg.add_to_pileup(
                                                k,
                                                PileupEntry::new(
                                                    PILEUP_ITEMS[if r { m + 4 } else { m }],
                                                    x_m,
                                                    x_h,
                                                ),
                                            )
                                        }
                                    } else if let Some(k) = pileup_ix {
                                        cpg.add_to_pileup(
                                            k,
                                            PileupEntry::new(
                                                if r {
                                                    PileupCode::UncalledRev
                                                } else {
                                                    PileupCode::UncalledFwd
                                                },
                                                x_m,
                                                x_h,
                                            ),
                                        )
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

pub(super) fn process_read_thread(
    cfg: &Config,
    ix: usize,
    contigs: &[String],
    reference: &Reference,
    out_send: Sender<CountBlock>,
    block_recv: Receiver<ProcessBlock>,
    block_send: Sender<BRecBlock>,
) -> anyhow::Result<()> {
    debug!("process read thread {} starting up", ix);

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
    debug!("process read thread {} shutting down", ix);
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
    let (tid, start) = get_read_coords(recs[0].brec())?;
    assert!(start <= end_pos);
    let (start, end_pos) = (start.saturating_sub(2), end_pos.saturating_add(2));
    
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
        end_pos,
    );
    let mut cb = CountBlock::new(idx, tid, start, rf, cfg, 2);

    for rec in recs {
        process_record(cfg, rf, rec, mm_parse, &mut cb)?
    }

    Ok(cb)
}
