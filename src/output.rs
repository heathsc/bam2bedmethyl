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
use crossbeam_channel::{Receiver, Sender, bounded};

use super::{
    config::Config,
    process_read::{
        count_block::CountBlock,
        cytosine::{CContext, Cytosine},
    },
};

use crate::read::pileup::PileupCode;

fn make_output_path(prefix: &str, context: &str) -> PathBuf {
    PathBuf::from(format!("{prefix}_{context}.bed"))
}

fn make_writer(prefix: &str, context: &str, ct: CompressType) -> anyhow::Result<Writer> {
    make_output_stream(make_output_path(prefix, context).as_path(), ct)
}

pub struct WriterBlock {
    cpg: Writer,
    chg: Option<Writer>,
    chh: Option<Writer>,
}

impl WriterBlock {
    fn new(prefix: &str, non_cpg: bool, ct: CompressType) -> anyhow::Result<Self> {
        let cpg = make_writer(prefix, "cpg", ct)?;
        let (chg, chh) = if non_cpg {
            (
                Some(make_writer(prefix, "chg", ct)?),
                Some(make_writer(prefix, "chh", ct)?),
            )
        } else {
            (None, None)
        };

        Ok(Self { cpg, chg, chh })
    }

    fn get_writer(&mut self, ctxt: CContext) -> Option<&mut Writer> {
        match ctxt {
            CContext::Cg => Some(&mut self.cpg),
            CContext::Chg(_) => self.chg.as_mut(),
            CContext::Chh(_,_) => self.chh.as_mut(),
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum WriterType {
    AllMeth,
    Mc,
    Hmc,
    HmFrac,
}

pub struct Writers {
    all_meth: Option<WriterBlock>,
    mc: Option<WriterBlock>,
    hmc: Option<WriterBlock>,
    hm_frac: Option<WriterBlock>,
    pileup: Option<Writer>,
}

impl Writers {
    fn new(cfg: &Config) -> anyhow::Result<Self> {
        let ct = if cfg.compress() {
            CompressType::Bgzip
        } else {
            CompressType::NoFilter
        };
        let non_cpg = cfg.non_cpg();
        let pref = cfg.output_prefix();

        let get_wb = |s| WriterBlock::new(format!("{pref}_{s}").as_str(), non_cpg, ct);

        let all_meth = Some(get_wb("meth")?);
        let mc = Some(get_wb("5mc")?);
        let hmc = Some(get_wb("5hmc")?);
        let hm_frac = Some(get_wb("hmc_frac")?);

        let pileup = if cfg.pileup() {
            Some(make_output_stream(
                &PathBuf::from(format!("{pref}_pipeline.tsv")),
                ct,
            )?)
        } else {
            None
        };

        Ok(Self {
            all_meth,
            mc,
            hmc,
            hm_frac,
            pileup,
        })
    }

    pub fn get_writer(&mut self, t: WriterType) -> WriterBlock {
        match t {
            WriterType::AllMeth => self.all_meth.take().unwrap(),
            WriterType::Mc => self.mc.take().unwrap(),
            WriterType::Hmc => self.hmc.take().unwrap(),
            WriterType::HmFrac => self.hm_frac.take().unwrap(),
        }
    }

    pub fn get_pipeline(&mut self) -> Option<Writer> {
        self.pileup.take()
    }
}

fn make_output_stream(path: &Path, ct: CompressType) -> anyhow::Result<Writer> {
    CompressIo::new()
        .path(path)
        .ctype(ct)
        .cthreads(CompressThreads::Set(4))
        .writer()
        .with_context(|| "Could not open output file/stream")
}

const OUTPUT_BLOCK_SIZE: usize = 8192;

struct OutputBlock {
    start: usize,
    tid: usize,
    cpgs: Vec<Cytosine>,
}

impl OutputBlock {
    fn init(start: usize, tid: usize) -> Self {
        Self {
            start,
            tid,
            cpgs: Vec::with_capacity(OUTPUT_BLOCK_SIZE),
        }
    }

    fn add_cpg(&mut self, mut cpg: Cytosine, delta: u32) -> bool {
        cpg.add_offset(delta);
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

    let wrt = Writers::new(cfg)?;

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
            trace!("Output thread received block {}", blk.idx());

            // Check if this is the block we are looking for
            if blk.idx() == idx {
                process_block(blk, &mut curr_blk, &mut out_blk, &block_send, ctg_names)?;
                idx += 1;
                // Check if we have the next blocks() in the pending hash
                while let Some(blk) = pending.remove(&idx) {
                    process_block(blk, &mut curr_blk, &mut out_blk, &block_send, ctg_names)?;
                    idx += 1;
                }
            } else {
                // Block has w=arrived out of order, so add to pending
                pending.insert(blk.idx(), blk);
            }
        }
        if let Some(cblk) = curr_blk.take() {
            flush_block_and_send(cblk, &mut out_blk, &block_send)?
        }
        drop(block_send);
        debug!("Output thread shutting down");
        jh.join().expect("Error joining output_block_thread")
    })
}

fn output_block_thread<'a>(
    s: &'a Scope<'a, '_>,
    mut wrt: Writers,
    ctg_names: &'a [String],
    r: Receiver<OutputBlock>,
) -> anyhow::Result<()> {
    debug!("Output block thread starting up");

    let mut jh = Vec::with_capacity(5);
    let mut chan = Vec::with_capacity(5);

    for wt in [
        WriterType::AllMeth,
        WriterType::Mc,
        WriterType::Hmc,
        WriterType::HmFrac,
    ] {
        let (snd, rcv) = bounded(2);
        chan.push(snd);
        let w = wrt.get_writer(wt);
        jh.push(s.spawn(move || writer_thread(w, ctg_names, rcv, wt)))
    }

    if let Some(w) = wrt.get_pipeline() {
        let (snd, rcv) = bounded(2);
        chan.push(snd);
        jh.push(s.spawn(move || pileup_writer_thread(w, ctg_names, rcv)))
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

fn writer_thread(
    mut w: WriterBlock,
    ctg_names: &[String],
    r: Receiver<Arc<OutputBlock>>,
    wt: WriterType,
) -> anyhow::Result<()> {
    debug!("Cpg writer thread {wt:?} starting up");

    match wt {
        WriterType::AllMeth => write_blocks(&mut w, ctg_names, "5mC+5hmC", r, |c| {
            [c.counts()[1], c.counts()[0]]
        })?,
        WriterType::Mc => write_blocks(&mut w, ctg_names, "5mC", r, |c| {
            [c.counts()[2], c.counts()[0] + c.counts()[3]]
        })?,
        WriterType::Hmc => write_blocks(&mut w, ctg_names, "5hmC", r, |c| {
            [c.counts()[3], c.counts()[0] + c.counts()[2]]
        })?,
        WriterType::HmFrac => write_blocks(&mut w, ctg_names, "5hmC/(5mC+5hmC)", r, |c| {
            [c.counts()[3], c.counts()[2]]
        })?,
    }
    debug!("Cpg writer thread {wt:?} shutting down");
    Ok(())
}

fn write_blocks<F: Fn(&Cytosine) -> [u32; 2]>(
    w: &mut WriterBlock,
    ctg_names: &[String],
    desc: &str,
    r: Receiver<Arc<OutputBlock>>,
    f: F,
) -> anyhow::Result<()> {
    for ob in r.iter() {
        let start = ob.start;
        let ctg = ctg_names[ob.tid].as_str();
        for c in ob.cpgs.iter() {
            write_cytosine(w, start, ctg, desc, c, &f)?
        }
    }
    Ok(())
}

fn write_cytosine<F: Fn(&Cytosine) -> [u32; 2]>(
    wb: &mut WriterBlock,
    start: usize,
    ctg: &str,
    desc: &str,
    c: &Cytosine,
    f: F,
) -> anyhow::Result<()> {
    let ct = f(c);
    if ct[0] + ct[1] > 0
        && let Some(w) = wb.get_writer(c.context())
    {
        let pos = (c.offset() as usize) + start;

        let cov = ct[0] + ct[1];
        let m = if cov > 0 {
            (ct[0] as f64) / (cov as f64)
        } else {
            0.0
        };

        let s1 = format!("{}\t{}", pos, pos + 1);
        writeln!(
            w,
            "{ctg}\t{s1}\t\"{desc}\"\t{}\t{}\t{s1}\t{}\t{cov}\t{:.4}\t{}",
            cov.min(1000),
            c.strand(),
            RGB_TAB[(m * 10.0 + 0.5) as usize],
            100.0 * m,
            c.context()
        )
        .with_context(|| "Error writing to bed file")
    } else {
        Ok(())
    }
}

fn pileup_writer_thread(
    mut w: Writer,
    ctg_names: &[String],
    r: Receiver<Arc<OutputBlock>>,
) -> anyhow::Result<()> {
    debug!("Pileup writer thread starting up");

    for ob in r.iter() {
        let start = ob.start;
        let ctg = ctg_names[ob.tid].as_str();
        for c in ob.cpgs.iter() {
            write_pileup(&mut w, start, ctg, c)?
        }
    }
    debug!("Pileup writer thread shutting down");
    Ok(())
}

fn write_pileup(w: &mut Writer, start: usize, ctg: &str, c: &Cytosine) -> anyhow::Result<()> {
    const ERR_STR: &str = "Error writing to pileup file";
    if let Some(v) = c.pileup() {
        let mut cts: [u32; 8] = [0; 8];
        for e in v.iter() {
            if let Some(ix) = match e.code() {
                PileupCode::CytosineFwd => Some(0),
                PileupCode::MethCytosineFwd => Some(1),
                PileupCode::HydroxyMethCytosineFwd => Some(2),
                PileupCode::TotalMethFwd => Some(3),
                PileupCode::CytosineRev => Some(4),
                PileupCode::MethCytosineRev => Some(5),
                PileupCode::HydroxyMethCytosineRev => Some(6),
                PileupCode::TotalMethRev => Some(7),
                _ => None,
            } {
                cts[ix] += 1
            }
        }

        if cts.iter().sum::<u32>() > 0 {
            let pos = (c.offset() as usize) + start + 1;
            write!(w, "{ctg}\t{pos}\t",).with_context(|| ERR_STR)?;
            for x in cts.iter() {
                write!(w, "{}\t", x).with_context(|| ERR_STR)?;
            }
            for e in v.iter() {
                write!(w, "{e}").with_context(|| ERR_STR)?;
            }
            writeln!(w).with_context(|| ERR_STR)?;
        }
    }
    Ok(())
}

fn process_block(
    mut blk: CountBlock,
    curr_blk: &mut Option<CountBlock>,
    out_blk: &mut Option<OutputBlock>,
    snd: &Sender<OutputBlock>,
    ctg_names: &[String],
) -> anyhow::Result<()> {
    trace!("Output thread processing block {}", blk.idx());

    if let Some(mut prev_blk) = curr_blk.take() {
        if prev_blk.tid() != blk.tid() {
            flush_block_and_send(prev_blk, out_blk, snd)?;
            debug!("Started output of {}", ctg_names[blk.tid()]);
        } else {
            assert!(blk.start() >= prev_blk.start());
            if blk.c_sites().is_empty() {
                blk = prev_blk
            } else {
                let delta = (blk.start() - prev_blk.start()) as u32;
                let first_x = blk.c_sites()[0].offset() + delta;
                write_partial_block(&mut prev_blk, out_blk, first_x, snd)?;
                blk.add_counts(&mut prev_blk);
            }
        }
    } else {
        debug!("Started output of {}", ctg_names[blk.tid()]);
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

fn get_ouput_block_and_delta(
    out_blk: &mut Option<OutputBlock>,
    blk: &CountBlock,
) -> (OutputBlock, u32) {
    if let Some(ob) = out_blk.take() {
        let delta = (blk.start() - ob.start) as u32;
        (ob, delta)
    } else {
        (OutputBlock::init(blk.start(), blk.tid()), 0)
    }
}

fn flush_block_and_send(
    mut blk: CountBlock,
    out_blk: &mut Option<OutputBlock>,
    snd: &Sender<OutputBlock>,
) -> anyhow::Result<()> {
    let (mut ob, delta) = get_ouput_block_and_delta(out_blk, &blk);

    for c in blk.c_sites_mut().drain(..) {
        ob.add_cpg(c, delta);
    }
    snd.send(ob).with_context(|| "Error sending output block")
}

fn write_partial_block(
    blk: &mut CountBlock,
    out_blk: &mut Option<OutputBlock>,
    first_x: u32,
    snd: &Sender<OutputBlock>,
) -> anyhow::Result<()> {
    let (mut ob, mut delta) = get_ouput_block_and_delta(out_blk, blk);

    let i = blk
        .c_sites()
        .binary_search_by_key(&first_x, |c| c.offset())
        .unwrap_or_else(|i| i);
    let (st, tid) = (blk.start(), blk.tid());

    for c in blk.c_sites_mut().drain(..i) {
        assert!(c.offset() < first_x);
        if ob.add_cpg(c, delta) {
            snd.send(ob).with_context(|| "Error sending output block")?;
            ob = OutputBlock::init(st, tid);
            delta = 0;
        }
    }
    *out_blk = Some(ob);
    Ok(())
}
