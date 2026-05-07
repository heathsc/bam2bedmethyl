use m_htslib::{hts::HtsPos, sam::BamRec};

use super::opt_index::OptIndex;

pub struct ReadRecord {
    brec: BamRec,
    end_pos: HtsPos,
    pileup_ix: OptIndex,
}

impl Default for ReadRecord {
    fn default() -> Self {
        Self {
            brec: BamRec::new(),
            end_pos: 0,
            pileup_ix: OptIndex::default(),
        }
    }
}

impl ReadRecord {
    #[inline]
    pub fn new() -> Self {
        Self::default()
    }

    #[inline]
    pub fn brec(&self) -> &BamRec {
        &self.brec
    }

    #[inline]
    pub fn brec_mut(&mut self) -> &mut BamRec {
        &mut self.brec
    }

    #[inline]
    pub fn end_pos(&self) -> HtsPos {
        self.end_pos
    }

    #[inline]
    pub fn pileup_ix(&self) -> Option<usize> {
        self.pileup_ix.get()
    }

    #[inline]
    pub(crate) fn update(&mut self, end_pos: HtsPos, pileup_ix: Option<usize>) {
        self.pileup_ix.set(pileup_ix);
        self.end_pos = end_pos;
    }
}
