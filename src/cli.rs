use std::{num::NonZeroUsize, path::PathBuf};

use clap::{command, value_parser, Arg, ArgAction, Command};

use super::utils::LogLevel;

pub fn cli_model() -> Command {
    command!()
        .arg(
            Arg::new("timestamp")
                .short('X')
                .long("timestamp")
                .value_parser(value_parser!(stderrlog::Timestamp))
                .value_name("GRANULARITY")
                .default_value("none")
                .help("Prepend log entries with a timestamp"),
        )
        .arg(
            Arg::new("loglevel")
                .short('l')
                .long("loglevel")
                .value_name("LOGLEVEL")
                .value_parser(value_parser!(LogLevel))
                .ignore_case(true)
                .default_value("info")
                .help("Set log level"),
        )
        .arg(
            Arg::new("quiet")
                .action(ArgAction::SetTrue)
                .long("quiet")
                .conflicts_with("loglevel")
                .help("Silence all output"),
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .value_parser(value_parser!(NonZeroUsize))
                .value_name("INT")
                .help("Set number of threads [default: available cores]"),
        )
        .arg(
            Arg::new("discard")
                .short('d')
                .long("discard-prob")
                .value_parser(value_parser!(f64))
                .value_name("PROB (0-1)")
                .default_value("0.0")
                .help("Prob. that reads will be (randomly) discarded for down sampling"),
        )
        .arg(
            Arg::new("seed")
                .short('S')
                .long("seed")
                .value_parser(value_parser!(u64))
                .value_name("SEED")
                .requires("discard")
                .help("Random number seed for down sampling"),
        )
        .arg(
            Arg::new("mapq_threshold")
                .short('q')
                .long("mapq-threshold")
                .value_parser(value_parser!(u8))
                .value_name("QUAL")
                .default_value("0")
                .help("Minimum mapping quality threshold"),
        )
        .arg(
            Arg::new("min_prob")
                .short('p')
                .long("min-prob")
                .value_parser(value_parser!(f64))
                .value_name("PROB (0.5-1)")
                .default_value("0.8")
                .help("Minimum prob. to call a methylation value"),
        )
        .arg(
            Arg::new("reference")
                .short('T')
                .long("reference")
                .value_parser(value_parser!(PathBuf))
                .value_name("REFERENCE_FILE")
                .required(true)
                .help("Input FASTA file with reference sequence"),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_parser(value_parser!(String))
                .value_name("OUTPUT_FILE")
                .default_value("bam2bedmethyl")
                .help("Prefix for output files"),
        )
        .arg(
            Arg::new("pileup")
                .action(ArgAction::SetTrue)
                .short('P')
                .long("pileup")
                .help("Generate methylation pileup)"),
        )
        .arg(
            Arg::new("non_cpg")
                .action(ArgAction::SetTrue)
                .short('n')
                .long("non-cpg")
                .help("Generate output for non-cpg as well as cpg contexts"),
        )
        .arg(
            Arg::new("compress")
                .action(ArgAction::SetTrue)
                .short('z')
                .long("compress")
                .help("Compress output (bgzip)"),
        )
        .arg(
            Arg::new("input")
                .value_parser(value_parser!(PathBuf))
                .value_name("INPUT_FILE")
                .help("Input file"),
        )
}
