use std::{
    num::NonZeroUsize,
    path::{Path, PathBuf},
};

use rand::{prelude::*, rngs::StdRng};

use rs_htslib::hts::HtsThreadPool;

use super::{cli::cli_model, utils::init_log};

#[derive(Debug)]
pub struct Config {
    // Input file - if none, input from stdin
    input_file: Option<PathBuf>,
    // Output file prefix - if not, set to bam2bedmethyl
    output_prefix: String,
    compress: bool,
    // Path to fasta file (possibly compressed) with reference sequence
    ref_file: PathBuf,
    // Down sampling - select proportion of reads that are discarded.  Must be in 0.0...=1.0
    discard: f64,
    seed: u64,
    // Filtering (carried out *after* any down sampling).  Note that prob threshold is scaled from 0
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
    pub fn output_prefix(&self) -> &str {
        &self.output_prefix
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

    let compress = m.get_flag("compress");
    let output_prefix = m
        .get_one::<String>("output")
        .map(|p| p.to_owned())
        .expect("Missing default output option");

    Ok(Config {
        input_file,
        output_prefix,
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
