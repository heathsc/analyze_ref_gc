use std::path::Path;

use anyhow::Context;
use compress_io::compress::CompressIo;
use serde::Serialize;

use crate::{betabin::write_hist, cli::Config, process::GcRes};

#[derive(Serialize)]
struct JsOutput<'a, 'b> {
    date: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    identifier: Option<&'a str>,
    #[serde(skip_serializing_if = "Option::is_none")]
    input: Option<&'a Path>,
    threads: usize,
    threshold: f64,
    read_lengths: &'a [u32],
    results: &'b GcRes,
}

impl<'a, 'b> JsOutput<'a, 'b> {
    fn make(cfg: &'a Config, results: &'b GcRes) -> Self {
        Self {
            date: cfg.date().to_rfc2822(),
            identifier: cfg.identifier(),
            input: cfg.input(),
            threads: cfg.threads(),
            threshold: cfg.threshold(),
            read_lengths: cfg.read_lengths(),
            results,
        }
    }
}

fn output_json<P: AsRef<Path>>(name: P, cfg: &Config, res: &GcRes) -> anyhow::Result<()> {
    debug!("Writing JSON output");
    let wrt = CompressIo::new()
        .path(name)
        .bufwriter()
        .with_context(|| "Could not open output JSON file")?;

    let out = JsOutput::make(cfg, res);

    serde_json::to_writer_pretty(wrt, &out)
        .with_context(|| "Error writing out JSON file with results")
}

fn output_dist<P: AsRef<Path>>(name: P, read_lengths: &[u32], res: &GcRes) -> anyhow::Result<()> {
    debug!("Writing expected GC distributions output");
    let mut wrt = CompressIo::new()
        .path(name)
        .bufwriter()
        .with_context(|| "Could not open output distribution file")?;

    write_hist(&mut wrt, read_lengths, res)
}

pub fn output(cfg: &Config, res: &GcRes) -> anyhow::Result<()> {
    let name = format!("{}.json", cfg.prefix());
    output_json(name, cfg, res)?;

    let name = format!("{}_dist.txt", cfg.prefix());
    output_dist(name, cfg.read_lengths(), res)
}
