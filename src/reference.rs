use std::{
    collections::{hash_map::Entry, HashMap},
    fmt::{self, Debug, Formatter},
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
    sync::Arc,
    thread,
};

use anyhow::Context;
use compress_io::compress::{CompressIo, Reader};
use crossbeam_channel::{unbounded, Receiver};
use rs_htslib::{
    faidx::{Faidx, Sequence},
    hts::{hts_get_log_level, hts_set_log_level, HtsLogLevel},
};

use crate::config::Config;

/// Load reference data from (possibly compressed) fasta file into Reference struct
pub fn load_reference(cfg: &Config) -> anyhow::Result<Reference> {
    let path = cfg.ref_file();
    info!("Loading reference from {:?}", path.display());

    // Open file
    let mut file = FastaFile::new(path)?;

    debug!("Reading in reference sequence from {}", path.display());

    let mut reference = Reference::new(path);
    while let Some(ctg) = file.read_fasta_record()? {
        reference.add_contig(ctg)?
    }

    info!("Finished loading reference");

    Ok(reference)
}

struct FastaFile<'a> {
    path: &'a Path,
    line: usize,
    buf: String,
    rdr: BufReader<Reader>,
}

impl<'a> FastaFile<'a> {
    fn new(path: &'a Path) -> anyhow::Result<Self> {
        let rdr = CompressIo::new()
            .path(path)
            .bufreader()
            .with_context(|| format!("Error opening FASTA file {}", path.display()))?;
        Ok(Self {
            path: path.as_ref(),
            line: 0,
            buf: String::new(),
            rdr,
        })
    }

    // Get next line from file.  Returns true if EOF
    fn get_line(&mut self) -> anyhow::Result<bool> {
        self.buf.clear();
        self.line += 1;
        Ok(self.rdr.read_line(&mut self.buf).with_context(|| {
            format!(
                "Error reading from {} at line {}",
                self.path.display(),
                self.line
            )
        })? == 0)
    }

    // Get next non-empty line.  If buffer non-empty returns buffer non-changed
    fn get_next_non_empty_line(&mut self) -> anyhow::Result<bool> {
        while self.buf.trim_end().is_empty() {
            if self.get_line()? {
                return Ok(true);
            }
        }
        Ok(false)
    }

    fn read_fasta_record(&mut self) -> anyhow::Result<Option<RefContig>> {
        loop {
            if self.get_next_non_empty_line()? {
                return Ok(None); // EOF
            }

            // Get contig name
            let name = Arc::from(parse_name(&self.buf).with_context(|| {
                format!(
                    "{}:{} Error reading FASTA name",
                    self.path.display(),
                    self.line
                )
            })?);

            // Read sequence
            trace!("Reading sequence {}", &name);
            let seq = self
                .read_sequence()
                .with_context(|| format!("Error reading sequence for contig {}", name))?;

            debug!("Read {}, length {}", &name, seq.len());
            return Ok(Some(RefContig {
                name,
                seq: RefSeq::Vec(seq),
            }));
        }
    }

    fn skip_sequence(&mut self) -> anyhow::Result<()> {
        while !self.get_line()? {
            if self.buf.starts_with('>') {
                break;
            }
            let s = self.buf.trim_end();
            if s.is_empty() {
                break;
            }
        }
        Ok(())
    }

    fn read_sequence(&mut self) -> anyhow::Result<Vec<u8>> {
        let mut seq = Vec::new();
        while !self.get_line()? {
            if self.buf.starts_with('>') {
                break;
            }
            let s = self.buf.trim_end();
            if s.is_empty() {
                break;
            }
            seq.extend_from_slice(s.as_bytes())
        }
        if seq.is_empty() {
            Err(anyhow!(
                "Empty sequence at {}:{}",
                self.path.display(),
                self.line
            ))
        } else {
            Ok(seq)
        }
    }
}

fn parse_name(s: &str) -> anyhow::Result<&str> {
    if let Some(t) = s.strip_prefix('>') {
        t.split_ascii_whitespace()
            .next()
            .ok_or_else(|| anyhow!("Missing name from FASTA record"))
    } else {
        Err(anyhow!(
            "Error parsing FASTA name: expecting '>', got '{}' ",
            s.chars().next().expect("Empty name string!")
        ))
    }
}

pub struct Reference {
    path: PathBuf,
    contig_seq: HashMap<Arc<str>, RefContig>,
}

impl Debug for Reference {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Reference [ path: {}, no: contigs: {} ]",
            self.path.display(),
            self.contig_seq.len()
        )
    }
}

impl Reference {
    fn new(p: &Path) -> Self {
        Self {
            path: p.to_owned(),
            contig_seq: HashMap::new(),
        }
    }

    pub fn get_seq(&self, ctg: &str, x: usize, y: usize) -> Option<&[u8]> {
        self.contig_seq.get(ctg).and_then(|c| c.seq(x, y))
    }

    fn add_contig(&mut self, ref_ctg: RefContig) -> anyhow::Result<()> {
        let ctg = ref_ctg.name.clone();
        match self.contig_seq.entry(ctg) {
            Entry::Occupied(_) => Err(anyhow!("Duplicate contig {} in fasta file", ref_ctg.name())),
            Entry::Vacant(e) => {
                e.insert(ref_ctg);
                Ok(())
            }
        }
    }
}

pub struct RefContig {
    name: Arc<str>,
    seq: RefSeq,
}

impl RefContig {
    pub fn name(&self) -> &str {
        self.name.as_ref()
    }

    pub fn seq(&self, x: usize, y: usize) -> Option<&[u8]> {
        self.seq.seq(x, y)
    }
}

pub enum RefSeq {
    Vec(Vec<u8>),
    Seq(Sequence),
}

impl RefSeq {
    pub fn seq(&self, x: usize, y: usize) -> Option<&[u8]> {
        match self {
            Self::Vec(v) => {
                assert!(x <= y);
                if x < v.len() {
                    let y = (y + 1).min(v.len());
                    Some(&v[x..y])
                } else {
                    None
                }
            }
            Self::Seq(s) => s.get_seq(x + 1, y + 1).ok(),
        }
    }
}

/// Try to load faidx index for reference file
fn try_load_faidx<P: AsRef<Path>>(path: P) -> Option<Faidx> {
    // Turn off hts error logs as we are not going to fail if we don't find the index
    let opt = hts_get_log_level();
    hts_set_log_level(HtsLogLevel::Off);

    let faidx = Faidx::load(path).ok();

    // Reset original hts logging level
    hts_set_log_level(opt);

    faidx
}

fn faidx_load_thread(
    ix: usize,
    cfg: &Config,
    mut f: Option<Faidx>,
    r: Receiver<usize>,
) -> anyhow::Result<Reference> {
    debug!("Reference read thread {} starting up", ix);
    let path = cfg.ref_file();
    let faidx = match f.take() {
        Some(x) => x,
        None => Faidx::load(path).with_context(|| {
            format!("Error loading index for reference file {}", path.display())
        })?,
    };
    let mut reference = Reference::new(path);
    for ctg_ix in r.iter() {
        let ctg = faidx.iseq(ctg_ix);

        debug!("({}) Reading in sequence from {}", ix, ctg);

        let seq = faidx
            .fetch_seq(&ctg, 0, None)
            .with_context(|| format!("Error reading sequencing for {}", ctg))?;

        let name = Arc::from(ctg.to_str()?);

        debug!("({}) Read sequence {}, length {}", ix, name, seq.len());

        reference.add_contig(RefContig {
            name,
            seq: RefSeq::Seq(seq),
        })?
    }
    debug!("Reference read thread {} shutting down", ix);
    Ok(reference)
}

/// Load reference from faidx indexed reference
fn load_reference_from_faidx(cfg: &Config, faidx: Faidx) -> anyhow::Result<Reference> {
    info!(
        "Loading faidx reference from {:?}",
        cfg.ref_file().display()
    );
    let nseq = faidx.nseq();
    let seq_ids: Vec<_> = (0..nseq).collect();
    let nt = cfg.threads().min(nseq);
    let mut faidx = Some(faidx);
    let (snd, r) = unbounded();
    let mut results = Vec::with_capacity(nt);
    thread::scope(|s| {
        let mut jh: Vec<_> = (0..nt)
            .map(|ix| {
                let f = faidx.take();
                let rc = r.clone();
                s.spawn(move || faidx_load_thread(ix, cfg, f, rc))
            })
            .collect();

        // Send contigs to process to readers
        for ctg in seq_ids {
            snd.send(ctg)
                .expect("Error sending new contig to SNV readers")
        }

        // Signal that no more contigs are coming so that the readers will exit
        drop(snd);

        // Wait for readers to finish
        for h in jh.drain(..) {
            let v = h.join().expect("Error while joining SNV reader");
            results.push(v);
        }
    });

    let mut reference = Reference::new(cfg.ref_file());
    for res in results.drain(..) {
        let mut rf = res?;
        for (_, v) in rf.contig_seq.drain() {
            reference.add_contig(v)?;
        }
    }
    info!("Finished loading reference");
    Ok(reference)
}

/// Check if reference file has faidx index.  If so then we read it in using htslib
/// (multithreaded if requested).  Otherwise we read it in as a standard (possibly compressed)
/// FASTA file.
pub fn handle_reference(cfg: &Config) -> anyhow::Result<Reference> {
    // Try to load faidx index
    if let Some(faidx) = try_load_faidx(cfg.ref_file()) {
        load_reference_from_faidx(cfg, faidx)
    } else {
        load_reference(cfg)
    }
}
