use rs_htslib::{hts::HtsPos, sam::BamRec};

use super::opt_index::OptIndex;

pub struct ReadRec {
    brec: BamRec,
    end_pos: HtsPos,
    pileup_ix: OptIndex,
}

impl Default for ReadRec {
    fn default() -> Self {
        Self {
            brec: BamRec::new().expect("Couldn't allocated new BamRec"),
            end_pos: 0,
            pileup_ix: OptIndex::default(),
        }
    }
}

impl ReadRec {
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
