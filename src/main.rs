#[macro_use]
extern crate anyhow;
#[macro_use]
extern crate log;

mod brec_block;
mod cli;
mod config;
mod output;
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
    let mut mm_parse = rs_htslib::sam::MMParse::default();
    read::read_input(&cfg, &r, &mut mm_parse)
}
