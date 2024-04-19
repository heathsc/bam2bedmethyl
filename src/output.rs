use anyhow::Context;
use std::{
    collections::HashMap,
    io::Write,
    path::{Path, PathBuf},
    sync::Arc,
    thread::{self, Scope},
};

use compress_io::{
    compress::{CompressIo, Writer},
    compress_type::{CompressThreads, CompressType},
};
use crossbeam_channel::{bounded, Receiver, Sender};

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
        .cthreads(CompressThreads::Set(2))
        .writer()
        .with_context(|| "Could not open output file/stream")
}

fn open_writers(cfg: &Config) -> anyhow::Result<[Option<Writer>; 3]> {
    let ct = if cfg.compress() {
        CompressType::Bgzip
    } else {
        CompressType::NoFilter
    };

    let output_paths = make_output_paths(cfg.output_prefix());

    Ok([
        Some(make_output_stream(&output_paths[0], ct)?),
        Some(make_output_stream(&output_paths[1], ct)?),
        Some(make_output_stream(&output_paths[2], ct)?),
    ])
}

const OUTPUT_BLOCK_SIZE: usize = 8192;

struct OutputBlock {
    start: usize,
    tid: usize,
    cpgs: Vec<CpG>,
}

impl OutputBlock {
    fn init(start: usize, tid: usize) -> Self {
        Self {
            start,
            tid,
            cpgs: Vec::with_capacity(OUTPUT_BLOCK_SIZE),
        }
    }

    fn add_cpg(&mut self, mut cpg: CpG, delta: u32) -> bool {
        cpg.offset += delta;
        self.cpgs.push(cpg);
        self.cpgs.len() >= OUTPUT_BLOCK_SIZE
    }
}

pub(super) fn output_thread(
    cfg: &Config,
    ctg_names: &[String],
    r: Receiver<CountBlock>,
) -> anyhow::Result<()> {
    debug!("Output thread starting up");

    let wrt = open_writers(cfg)?;

    // Get channels for output blocks
    let (block_send, block_recv) = bounded(6);

    thread::scope(|s| {
        // Spawn output block thread (Combined)
        let jh = s.spawn(|| output_block_thread(s, wrt, ctg_names, block_recv));

        // Index of next block for output
        let mut idx = 0;

        // Current block being output
        let mut curr_blk: Option<CountBlock> = None;
        let mut out_blk: Option<OutputBlock> = None;

        // Map containing blocks that arrived out of order
        let mut pending = HashMap::new();

        for blk in r.iter() {
            trace!("Output thread received block {}", blk.idx);

            // Check if this is the block we are looking for
            if blk.idx == idx {
                process_block(blk, &mut curr_blk, &mut out_blk, &block_send, ctg_names)?;
                idx += 1;
                // Check if we have the next blocks() in the pending hash
                while let Some(blk) = pending.remove(&idx) {
                    process_block(blk, &mut curr_blk, &mut out_blk, &block_send, ctg_names)?;
                    idx += 1;
                }
            } else {
                // Block has w=arrived out of order, so add to pending
                pending.insert(blk.idx, blk);
            }
        }
        if let Some(cblk) = curr_blk.take() {
            flush_block_and_send(&cblk, &mut out_blk, &block_send)?
        }
        drop(block_send);
        debug!("Output thread shutting down");
        jh.join().expect("Error joining output_block_thread")
    })
}

fn output_block_thread<'a>(
    s: &'a Scope<'a, '_>,
    mut wrt: [Option<Writer>; 3],
    ctg_names: &'a [String],
    r: Receiver<OutputBlock>,
) -> anyhow::Result<()> {
    debug!("Output block thread starting up");

    let mut jh = Vec::with_capacity(3);
    let mut chan = Vec::with_capacity(3);

    for (ix, w1) in wrt.iter_mut().enumerate() {
        let (snd, rcv) = bounded(2);
        chan.push(snd);
        let w = w1.take().unwrap();
        jh.push(s.spawn(move || cpg_writer_thread(w, ctg_names, rcv, ix)))
    }

    for ob in r.iter() {
        let ob = Arc::new(ob);
        for c in chan.iter() {
            let ob_clone = ob.clone();
            c.send(ob_clone)
                .expect("Error sending block to cpg writer threads")
        }
    }

    drop(chan);
    debug!("Output block thread shutting down");
    for j in jh.drain(..) {
        j.join().expect("Error joining cpg writer threads")?
    }
    Ok(())
}

fn cpg_writer_thread(
    w: Writer,
    ctg_names: &[String],
    r: Receiver<Arc<OutputBlock>>,
    ix: usize,
) -> anyhow::Result<()> {
    debug!("Cpg writer thread {ix} starting up");

    match ix {
        0 => write_cpg_blocks(w, ctg_names, r, |c| {
            (
                [c.fwd_counts[1], c.fwd_counts[0]],
                [c.rev_counts[1], c.rev_counts[0]],
                "5mC+5hmC",
            )
        })?,
        1 => write_cpg_blocks(w, ctg_names, r, |c| {
            (
                [c.fwd_counts[2], c.fwd_counts[0] + c.fwd_counts[3]],
                [c.rev_counts[2], c.rev_counts[0] + c.rev_counts[3]],
                "5mC",
            )
        })?,
        2 => write_cpg_blocks(w, ctg_names, r, |c| {
            (
                [c.fwd_counts[3], c.fwd_counts[0] + c.fwd_counts[2]],
                [c.rev_counts[3], c.rev_counts[0] + c.rev_counts[2]],
                "5hmC",
            )
        })?,
        _ => panic!("Illegal options"),
    }
    debug!("Cpg writer thread {ix} shutting down");
    Ok(())
}

fn write_cpg_blocks<F: Fn(&CpG) -> ([u32; 2], [u32; 2], &str)>(
    mut w: Writer,
    ctg_names: &[String],
    r: Receiver<Arc<OutputBlock>>,
    f: F,
) -> anyhow::Result<()> {
    for ob in r.iter() {
        let start = ob.start;
        let ctg = ctg_names[ob.tid].as_str();
        for c in ob.cpgs.iter() {
            write_cpg(&mut w, start, ctg, c, &f)?
        }
    }
    Ok(())
}

fn write_cpg<F: Fn(&CpG) -> ([u32; 2], [u32; 2], &str)>(
    w: &mut Writer,
    start: usize,
    ctg: &str,
    c: &CpG,
    f: F,
) -> anyhow::Result<()> {
    let (forward, reverse, desc) = f(c);
    if forward[0] + forward[1] + reverse[0] + reverse[1] > 0 {
        let pos = (c.offset as usize) + start;

        let fm = |ct: [u32; 2]| -> (u32, f64) {
            let cov = ct[0] + ct[1];
            let m = if cov > 0 {
                (ct[0] as f64) / (cov as f64)
            } else {
                0.0
            };
            (cov, m)
        };

        let (cov, m) = fm(forward);
        let s1 = format!("{}\t{}", pos, pos + 1);
        writeln!(
            w,
            "{ctg}\t{s1}\t\"{desc}\"\t{}\t-\t{s1}\t{}\t{cov}\t{:.2}",
            cov.min(1000),
            RGB_TAB[(m * 10.0 + 0.5) as usize],
            100.0 * m
        )
        .with_context(|| "Error writing to bed file")?;
        let (cov, m) = fm(reverse);
        let s1 = format!("{}\t{}", pos + 1, pos + 2);
        writeln!(
            w,
            "{ctg}\t{s1}\t\"{desc}\"\t{}\t-\t{s1}\t{}\t{cov}\t{:.2}",
            cov.min(1000),
            RGB_TAB[(m * 10.0 + 0.5) as usize],
            100.0 * m
        )
        .with_context(|| "Error writing to bed file")
    } else {
        Ok(())
    }
}

fn process_block(
    mut blk: CountBlock,
    curr_blk: &mut Option<CountBlock>,
    out_blk: &mut Option<OutputBlock>,
    snd: &Sender<OutputBlock>,
    ctg_names: &[String],
) -> anyhow::Result<()> {
    trace!("Output thread processing block {}", blk.idx);

    if let Some(mut prev_blk) = curr_blk.take() {
        if prev_blk.tid != blk.tid {
            flush_block_and_send(&prev_blk, out_blk, snd)?;
            debug!("Started output of {}", ctg_names[blk.tid]);
        } else {
            assert!(blk.start >= prev_blk.start);
            if blk.cpg_sites.is_empty() {
                blk = prev_blk
            } else {
                let delta = (blk.start - prev_blk.start) as u32;
                let first_x = blk.cpg_sites[0].offset + delta;
                write_partial_block(&prev_blk, out_blk, first_x, snd)?;
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

fn flush_block_and_send(
    blk: &CountBlock,
    out_blk: &mut Option<OutputBlock>,
    snd: &Sender<OutputBlock>,
) -> anyhow::Result<()> {
    let (mut ob, delta) = if let Some(ob) = out_blk.take() {
        let delta = (blk.start - ob.start) as u32;
        (ob, delta)
    } else {
        (OutputBlock::init(blk.start, blk.tid), 0)
    };
    for c in blk.cpg_sites.iter() {
        ob.add_cpg(*c, delta);
    }
    snd.send(ob).with_context(|| "Error sending output block")
}
fn write_partial_block(
    blk: &CountBlock,
    out_blk: &mut Option<OutputBlock>,
    first_x: u32,
    snd: &Sender<OutputBlock>,
) -> anyhow::Result<()> {
    let (mut ob, mut delta) = if let Some(ob) = out_blk.take() {
        let delta = (blk.start - ob.start) as u32;
        (ob, delta)
    } else {
        (OutputBlock::init(blk.start, blk.tid), 0)
    };
    for c in blk.cpg_sites.iter() {
        if c.offset >= first_x {
            break;
        }
        if ob.add_cpg(*c, delta) {
            snd.send(ob).with_context(|| "Error sending output block")?;
            ob = OutputBlock::init(blk.start, blk.tid);
            delta = 0;
        }
    }
    *out_blk = Some(ob);
    Ok(())
}
