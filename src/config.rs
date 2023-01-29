use std::{
    num::NonZeroUsize,
    path::{Path, PathBuf},
};

use clap::ArgMatches;
use rand::{prelude::*, rngs::StdRng};

use crate::utils;
use rs_htslib::hts::HtsThreadPool;

use super::{cli::cli_model, utils::init_log};

#[derive(Debug)]
pub struct Config {
    // Input file - if none, input from stdin
    input_file: Option<PathBuf>,
    // Output file - if not, output to stdout
    output_file: Option<PathBuf>,
    compress: bool,
    // Path to fasta file (possibly compressed) with reference sequence
    ref_file: PathBuf,
    // Down sampling - select proportion of reads that are discarded.  Must be in 0.0..=1.0
    discard: f64,
    seed: u64,
    // Filtering (carried out *after* any down sampling.  Note that prob threshold is scaled from 0
    // to 255 so that x is a prob. between x / 256 and (x + 1) / 256
    mapq_threshold: u8,
    prob_threshold: u8,
    // General options:
    threads: usize,
    hts_thread_pool: Option<HtsThreadPool>,
}

impl Config {
    pub fn input_file(&self) -> Option<&Path> {
        self.input_file.as_deref()
    }
    pub fn output_file(&self) -> Option<&Path> {
        self.output_file.as_deref()
    }
    pub fn compress(&self) -> bool {
        self.compress
    }
    pub fn ref_file(&self) -> &Path {
        &self.ref_file.as_path()
    }
    pub fn discard(&self) -> f64 {
        self.discard
    }
    pub fn seed(&self) -> u64 {
        self.seed
    }
    pub fn mapq_threshold(&self) -> u8 {
        self.mapq_threshold
    }
    pub fn prob_threshold(&self) -> u8 {
        self.prob_threshold
    }
    pub fn threads(&self) -> usize {
        self.threads
    }
    pub fn hts_thread_pool(&self) -> Option<&HtsThreadPool> {
        self.hts_thread_pool.as_ref()
    }
}

/// Get output file and compress options from command line
/// If compress is not set and the output file is set, set compress option based on the output
/// filename extension.  
/// If compress is set and output filename is set, add a gz extension if not already present.
/// If compress is set and output filename is not set, check whether we are outputting to a tty
/// and, if so, turn off compression
fn get_output(m: &ArgMatches) -> (Option<PathBuf>, bool) {
    // Get options from CLI
    let mut compress = m.get_flag("compress");
    let mut pathbuf = m.get_one::<PathBuf>("output").map(|p| p.to_owned());

    if let Some(mut p) = pathbuf.take() {
        // If output filename set...
        if let Some(ext) = p.extension() {
            // If filename has an extension...
            if ext == "gz" {
                // If extension is gz, set compress irrespective of CLI compress option
                compress = true
            } else if compress {
                // Otherwise if CLI compress option set, add gz extension to filename
                let mut s = p.into_os_string();
                s.push(".gz");
                p = PathBuf::from(s)
            }
        } else if compress {
            // If CLI option set and the filename has no extension, add gz extension to filename
            p.set_extension("gz");
        }
        pathbuf = Some(p)
    } else {
        // No output filename set.  Check if compress option is set and output is to tty.
        // If so, turn off compression
        if compress && utils::isatty(libc::STDOUT_FILENO) {
            warn!("Terminal output will not be compressed");
            compress = false;
        }
    }

    (pathbuf, compress)
}

pub fn handle_cli() -> anyhow::Result<Config> {
    // Get matches from command line
    let m = cli_model().get_matches();

    // Setup logging

    init_log(&m);

    debug!("Processing command line options");

    let ref_file = m
        .get_one::<PathBuf>("reference")
        .expect("Missing reference")
        .to_owned();

    let input_file = m.get_one::<PathBuf>("input").map(|p| p.to_owned());

    let discard = match m
        .get_one::<f64>("discard")
        .expect("Missing default discard value")
    {
        x if (0.0..=1.0).contains(x) => Ok(*x),
        x => Err(anyhow!("discard option {} not in range 0-1", x)),
    }?;

    // If seed is not given then we generate a random seed, using the default RNG seeded from system genrand()
    let seed = m.get_one::<u64>("seed").copied().unwrap_or_else(|| {
        let mut rng = StdRng::from_entropy();
        rng.next_u64()
    });

    let mapq_threshold = *m
        .get_one::<u8>("mapq_threshold")
        .expect("Missing default value");

    let prob_threshold = match m.get_one::<f64>("min_prob").expect("Missing default value") {
        x if (0.5..=1.0).contains(x) => Ok((x * 256.0).round().min(255.0) as u8),
        x => Err(anyhow!("meth probability option {} not in range 0.5-1", x)),
    }?;

    // Threads option should be non-zero.  If not set, set to number of available CPUs
    let threads = m
        .get_one::<NonZeroUsize>("threads")
        .map(|i| usize::from(*i))
        .unwrap_or_else(|| num_cpus::get());

    // If multiple treads requested, set up HtsThreadPool
    let hts_thread_pool = if threads > 1 {
        HtsThreadPool::new(threads)
    } else {
        None
    };

    let (output_file, compress) = get_output(&m);

    Ok(Config {
        input_file,
        output_file,
        ref_file,
        compress,
        discard,
        seed,
        mapq_threshold,
        prob_threshold,
        threads,
        hts_thread_pool,
    })
}
