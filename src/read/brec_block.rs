use rs_htslib::{hts::HtsRead, sam::SamReader};

use super::read_record::ReadRec;

/// Storage for consecutive BAM records read from
/// input file/stream
///
/// brec_vec is always fully populated with BamRec, and ix is the index of the next available
/// element.  If ix == brec_vec.len() then the block is full and nothing can be added
const BREC_BLOCK_SIZE: usize = 256;

pub struct BRecBlock {
    ix: usize,
    read_rec_vec: Vec<ReadRec>,
}

impl Default for BRecBlock {
    fn default() -> Self {
        let mut brec_vec = Vec::with_capacity(BREC_BLOCK_SIZE);
        for _ in 0..BREC_BLOCK_SIZE {
            brec_vec.push(ReadRec::new())
        }
        Self {
            ix: 0,
            read_rec_vec: brec_vec,
        }
    }
}

impl BRecBlock {
    #[inline]
    pub fn new() -> Self {
        Self::default()
    }

    #[inline]
    pub(super) fn clear(&mut self) {
        self.ix = 0
    }

    pub(super) fn next_rec<'a>(
        &'a mut self,
        rdr: &mut SamReader,
        pending: &'a mut Option<ReadRec>,
    ) -> anyhow::Result<BrecFill<'a>> {
        Ok(match self.read_rec_vec.get_mut(self.ix) {
            Some(b) => {
                if let Some(r) = pending.take() {
                    assert_eq!(self.ix, 0, "Pending on non-empty block");
                    let _ = std::mem::replace(b, r);
                    self.ix += 1;
                    BrecFill::Rec(b)
                } else if !rdr.read(b.brec_mut())? {
                    BrecFill::EndOfFile
                } else {
                    self.ix += 1;
                    BrecFill::Rec(b)
                }
            }
            None => BrecFill::EndOfBlock,
        })
    }

    #[inline]
    pub(super) fn push_back(&mut self) -> ReadRec {
        self.decr_ix();
        let mut r = ReadRec::new();
        std::mem::swap(&mut r, &mut self.read_rec_vec[self.ix]);
        r
    }

    #[inline]
    pub(super) fn decr_ix(&mut self) {
        assert!(self.ix > 0);
        self.ix -= 1;
    }

    #[inline]
    pub fn brec_vec(&mut self) -> &mut [ReadRec] {
        &mut self.read_rec_vec[..self.ix]
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.ix == 0
    }
}

pub(super) enum BrecFill<'a> {
    Rec(&'a mut ReadRec),
    EndOfBlock,
    EndOfFile,
}
