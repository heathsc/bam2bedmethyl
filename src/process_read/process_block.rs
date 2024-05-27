use crate::read::brec_block::BRecBlock;

pub(crate) struct ProcessBlock {
    pub(super) idx: usize,
    pub(super) end_pos: usize,
    pub(super) bblock: BRecBlock,
}

impl ProcessBlock {
    pub(crate) fn new(bblock: BRecBlock, idx: usize, end_pos: usize) -> Self {
        Self {
            idx,
            end_pos,
            bblock,
        }
    }
}
