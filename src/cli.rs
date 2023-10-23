use std::path::{Path, PathBuf};

use chrono::{DateTime, Local};

mod cli_model;

pub struct Config {
    input: Option<PathBuf>,
    prefix: String,
    identifier: Option<String>,
    threads: usize,
    threshold: f64,
    bisulfite: bool,
    read_lengths: Vec<u32>,
    date: DateTime<Local>,
}

impl Config {
    pub fn input(&self) -> Option<&Path> {
        self.input.as_deref()
    }

    pub fn threads(&self) -> usize {
        self.threads
    }

    pub fn read_lengths(&self) -> &[u32] {
        &self.read_lengths
    }

    pub fn threshold(&self) -> f64 {
        self.threshold
    }

    pub fn prefix(&self) -> &str {
        self.prefix.as_str()
    }

    pub fn identifier(&self) -> Option<&str> {
        self.identifier.as_deref()
    }
    
    pub fn date(&self) -> &DateTime<Local> { &self.date }
    
    pub fn bisulfite(&self) -> bool { self.bisulfite }
}

pub fn handle_cli() -> anyhow::Result<Config> {
    let c = cli_model::cli_model();
    let m = c.get_matches();
    super::utils::init_log(&m);

    let input = m.get_one::<PathBuf>("input").map(|p| p.to_owned());

    let threads = m
        .get_one::<u64>("threads")
        .map(|x| *x as usize)
        .unwrap_or_else(num_cpus::get);

    let read_lengths: Vec<u32> = m
        .get_many("read_lengths")
        .expect("Missing default argument")
        .copied()
        .collect();

    let threshold = match m
        .get_one::<f64>("threshold")
        .expect("Missing default argument")
    {
        x if x > &0.0 && x <= &1.0 => Ok(*x),
        _ => Err(anyhow!("Illegal threshold: must be > 0 and <= 1.0")),
    }?;

    let prefix = m
        .get_one::<String>("prefix")
        .map(|s| s.to_owned())
        .expect("Missing default argument");

    let identifier = m.get_one::<String>("identifier").map(|s| s.to_owned());

    let bisulfite = !m.get_flag("no_bisulfite");
    
    Ok(Config {
        input,
        prefix,
        identifier,
        threads,
        bisulfite,
        threshold,
        read_lengths,
        date: Local::now(),
    })
}
