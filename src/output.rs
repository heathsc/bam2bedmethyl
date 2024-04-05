use anyhow::Context;
use std::{
    collections::HashMap,
    io::Write,
    path::{Path, PathBuf},
};

use compress_io::{
    compress::{CompressIo, Writer},
    compress_type::{CompressThreads, CompressType},
};
use crossbeam_channel::Receiver;

use super::{
    config::Config,
    read::{CountBlock, CpG},
};

fn make_output_paths(prefix: &str) -> [PathBuf; 3] {
    let path_combined = PathBuf::from(format!("{prefix}_cpg.bed"));
    let path_meth = PathBuf::from(format!("{prefix}_5mC_cpg.bed"));
    let path_hmeth = PathBuf::from(format!("{prefix}_5hmC_cpg.bed"));
    [path_combined, path_meth, path_hmeth]
}

fn make_output_stream(path: &Path, ct: CompressType) -> anyhow::Result<Writer> {
    CompressIo::new()
        .path(path)
        .ctype(ct)
        .cthreads(CompressThreads::NPhysCores)
        .writer()
        .with_context(|| "Could not open output file/stream")
}

fn open_writers(cfg: &Config) -> anyhow::Result<[Writer; 3]> {
    let ct = if cfg.compress() {
        CompressType::Bgzip
    } else {
        CompressType::NoFilter
    };

    let output_paths = make_output_paths(cfg.output_prefix());

    Ok([
        make_output_stream(&output_paths[0], ct)?,
        make_output_stream(&output_paths[1], ct)?,
        make_output_stream(&output_paths[2], ct)?,
    ])
}

pub(super) fn output_thread(
    cfg: &Config,
    ctg_names: &[String],
    r: Receiver<CountBlock>,
) -> anyhow::Result<()> {
    debug!("Output thread starting up");

    let mut wrt = open_writers(cfg)?;

    // Index of next block for output
    let mut idx = 0;

    // Current block being output
    let mut curr_blk: Option<CountBlock> = None;

    // Map containing blocks that arrived out of order
    let mut pending = HashMap::new();

    for blk in r.iter() {
        trace!("Output thread received block {}", blk.idx);

        // Check if this is the block we are looking for
        if blk.idx == idx {
            output_block(&mut wrt, blk, &mut curr_blk, ctg_names)?;
            idx += 1;
            // Check if we have the next blocks() in the pending hash
            idx = output_pending(&mut wrt, &mut pending, idx, &mut curr_blk, ctg_names)?;
        } else {
            // Block has w=arrived out of order, so add to pending
            pending.insert(blk.idx, blk);
        }
    }

    debug!("Output thread shutting down");
    Ok(())
}

fn output_block(
    wrt: &mut [Writer; 3],
    mut blk: CountBlock,
    curr_blk: &mut Option<CountBlock>,
    ctg_names: &[String],
) -> anyhow::Result<()> {
    trace!("Output thread outputting block {}", blk.idx);

    if let Some(mut prev_blk) = curr_blk.take() {
        if prev_blk.tid != blk.tid {
            flush_block(wrt, &prev_blk, ctg_names[prev_blk.tid].as_str())?;
            debug!("Started output of {}", ctg_names[blk.tid]);
        } else {
            assert!(blk.start >= prev_blk.start);
            if blk.cpg_sites.is_empty() {
                blk = prev_blk
            } else {
                let delta = (blk.start - prev_blk.start) as u32;
                let first_x = blk.cpg_sites[0].offset + delta;
                write_partial_block(wrt, &prev_blk, first_x, ctg_names[prev_blk.tid].as_str())?;
                blk.add_counts(&mut prev_blk);
            }
        }
    } else {
        debug!("Started output of {}", ctg_names[blk.tid]);
    }
    *curr_blk = Some(blk);
    Ok(())
}

const RGB_TAB: [&str; 11] = [
    "0,255,0",
    "55,255,0",
    "105,255,0",
    "155,255,0",
    "205,255,0",
    "255,255,0",
    "255,205,0",
    "255,155,0",
    "255,105,0",
    "255,55,0",
    "255,0,0",
];

fn write_bed_entry(
    wrt: &mut Writer,
    ctg: &str,
    pos: usize,
    strand: char,
    ct_mod: u32,
    ct_nmod: u32,
) -> anyhow::Result<()> {
    let cov = ct_nmod + ct_mod;
    let m = if cov > 0 {
        (ct_mod as f64) / (cov as f64)
    } else {
        0.0
    };
    writeln!(
        wrt,
        "{}\t{}\t{}\t\".\"\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}",
        ctg,
        pos,
        pos + 1,
        cov.min(1000),
        strand,
        pos,
        pos + 1,
        RGB_TAB[(m * 10.0 + 0.5) as usize],
        cov,
        100.0 * m
    )
    .with_context(|| "Error writing output")
}

fn write_bed_line(
    wrt: &mut [Writer; 3],
    ctg: &str,
    pos: usize,
    strand: char,
    cts: &[u32],
) -> anyhow::Result<()> {
    write_bed_entry(&mut wrt[0], ctg, pos, strand, cts[1], cts[0])?;
    write_bed_entry(&mut wrt[1], ctg, pos, strand, cts[2], cts[3] + cts[0])?;
    write_bed_entry(&mut wrt[2], ctg, pos, strand, cts[3], cts[2] + cts[0])
}

fn write_cpg_entry(wrt: &mut [Writer; 3], c: &CpG, start: usize, ctg: &str) -> anyhow::Result<()> {
    if c.fwd_counts[0] + c.fwd_counts[1] + c.rev_counts[0] + c.rev_counts[1] > 0 {
        let pos = (c.offset as usize) + start;
        write_bed_line(wrt, ctg, pos, '+', &c.fwd_counts)?;
        write_bed_line(wrt, ctg, pos + 1, '-', &c.rev_counts)?;
    }
    Ok(())
}

fn flush_block(wrt: &mut [Writer; 3], blk: &CountBlock, ctg: &str) -> anyhow::Result<()> {
    for c in blk.cpg_sites.iter() {
        write_cpg_entry(wrt, c, blk.start, ctg)?
    }
    Ok(())
}

fn write_partial_block(
    wrt: &mut [Writer; 3],
    blk: &CountBlock,
    first_x: u32,
    ctg: &str,
) -> anyhow::Result<()> {
    for c in blk.cpg_sites.iter() {
        if c.offset >= first_x {
            break;
        }
        write_cpg_entry(wrt, c, blk.start, ctg)?
    }
    Ok(())
}

fn output_pending(
    wrt: &mut [Writer; 3],
    pending: &mut HashMap<usize, CountBlock>,
    mut idx: usize,
    curr_blk: &mut Option<CountBlock>,
    ctg_names: &[String],
) -> anyhow::Result<usize> {
    while let Some(blk) = pending.remove(&idx) {
        output_block(wrt, blk, curr_blk, ctg_names)?;
        idx += 1;
    }
    Ok(idx)
}
