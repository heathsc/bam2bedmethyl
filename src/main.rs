#[macro_use]
extern crate anyhow;
#[macro_use]
extern crate log;

mod cli;
mod config;
mod count_block;
mod output;
mod process_read;
mod read;
mod reference;
mod utils;

fn main() -> anyhow::Result<()> {
    // Set up configuration from CLU
    let cfg = config::handle_cli()?;
    debug!("{:?}", cfg);

    // Read in reference file
    let r = reference::handle_reference(&cfg)?;

    // Process input
    read::read_input(&cfg, &r)
}
