use std::path::PathBuf;

use clap::{command, value_parser, Arg, ArgAction, Command};

use crate::utils::LogLevel;

pub(super) fn cli_model() -> Command {
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
                .value_parser(value_parser!(u64).range(1..))
                .value_name("INT")
                .help("Set number of process threads [default: number of available cores]"),
        )
        .arg(
            Arg::new("threshold")
                .short('T')
                .long("threshold")
                .value_parser(value_parser!(f64))
                .value_name("PROPORTION")
                .default_value("0.8")
                .help("Set threshold (0 > x <= 1) for proportion of bases required"),
        )
        .arg(
            Arg::new("prefix")
                .short('p')
                .long("prefix")
                .value_parser(value_parser!(String))
                .value_name("PREFIX")
                .default_value("analyze_gc")
                .help("Set prefix for output file names"),
        )
        .arg(
            Arg::new("identifier")
                .short('i')
                .long("identifier")
                .value_parser(value_parser!(String))
                .value_name("IDENTIFIER")
                .help("Set identifier for reference genome"),
        )
        .arg(
            Arg::new("read_lengths")
                .short('r')
                .long("read_lengths")
                .value_parser(value_parser!(u32).range(1..))
                .value_name("INT")
                .num_args(1..)
                .default_values(["50", "75", "100", "150", "200", "250", "300"])
                .help("Set read lengths to analyze"),
        )
        .arg(
            Arg::new("input")
                .value_parser(value_parser!(PathBuf))
                .value_name("INPUT")
                .help("Input FASTA file"),
        )
}
