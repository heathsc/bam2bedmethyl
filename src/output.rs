use std::collections::HashMap;

use crossbeam_channel::{Receiver, Sender};

use super::{brec_block::BRecBlock, config::Config, read::CountBlock};

pub(super) fn output_thread(
    cfg: &Config,
    ctg_names: &[String],
    r: Receiver<CountBlock>,
) -> anyhow::Result<()> {
    debug!("Output thread starting up");
    for blk in r.iter() {}

    debug!("Output thread shutting down");
    Ok(())
}
