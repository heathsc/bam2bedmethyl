use rs_htslib::sam::BamRec;

/// Storage for consecutive BAM records read from
/// input file/stream
///
/// brec_vec is always fully populated with BamRec, and ix is the index of the next available
/// element.  If ix == brec_vec.len() then the block is full and nothing can be added
///

const BREC_BLOCK_SIZE: usize = 256;

pub(super) struct BRecBlock {
    ix: usize,
    brec_vec: Vec<BamRec>,
}

impl Default for BRecBlock {
    fn default() -> Self {
        let mut brec_vec = Vec::with_capacity(BREC_BLOCK_SIZE);
        for _ in 0..BREC_BLOCK_SIZE {
            brec_vec.push(BamRec::new().expect("Couldn't allocate Bam Records"))
        }
        Self { ix: 0, brec_vec }
    }
}

impl BRecBlock {
    pub(super) fn new() -> Self {
        Self::default()
    }

    pub(super) fn clear(&mut self) {
        self.ix = 0
    }

    pub(super) fn next_rec(&mut self) -> Option<&mut BamRec> {
        match self.brec_vec.get_mut(self.ix) {
            Some(b) => {
                self.ix += 1;
                Some(b)
            }
            None => None,
        }
    }

    pub(super) fn decr_ix(&mut self) {
        assert!(self.ix > 0);
        self.ix -= 1;
    }

    pub(super) fn brec_vec(&mut self) -> &mut [BamRec] {
        &mut self.brec_vec[..self.ix]
    }

    pub(super) fn is_empty(&self) -> bool {
        self.ix == 0
    }
}

pub(super) struct ProcessBlock {
    pub(super) idx: usize,
    pub(super) end_pos: usize,
    pub(super) bblock: BRecBlock,
}

impl ProcessBlock {
    pub(super) fn new(bblock: BRecBlock, idx: usize, end_pos: usize) -> Self {
        Self {
            idx,
            end_pos,
            bblock,
        }
    }
}
