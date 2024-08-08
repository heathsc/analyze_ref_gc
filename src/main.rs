#[macro_use]
extern crate log;
#[macro_use]
extern crate anyhow;

mod betabin;
mod cli;
mod output;
mod process;
mod reader;
mod utils;
mod regions;

fn main() -> anyhow::Result<()> {
    let cfg = cli::handle_cli()?;
    let res = process::process(&cfg)?;
    output::output(&cfg, &res)
}
