use std::fmt;

use crate::read::pileup::PileupEntry;

#[derive(Debug, Copy, Clone)]
pub enum CContext {
    Cg,
    Chg(u8),
    Chh(u8, u8),
}

impl fmt::Display for CContext {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Cg => write!(f, "CG"),
            Self::Chg(c) => write!(f, "C{}G", *c as char),
            Self::Chh(c1,c2) => write!(f, "C{}{}", *c1 as char, *c2 as char),
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Forward => "+",
                Self::Reverse => "-",
            }
        )
    }
}

/// The counts are:
///   0 - not modified
///   1 - 5mC or 5hmC
///   2 - 5mC
///   3 - 5hmC
#[derive(Debug)]
pub(crate) struct Cytosine {
    offset: u32,
    strand: Strand,
    context: CContext,
    counts: [u32; 4],
    pileup: Option<Vec<PileupEntry>>,
}

impl Cytosine {
    #[inline]
    pub(super) fn new(offset: u32, strand: Strand, context: CContext, gen_pileup: bool) -> Self {
        Self {
            offset,
            strand,
            context,
            counts: [0; 4],
            pileup: if gen_pileup { Some(Vec::new()) } else { None },
        }
    }

    #[inline]
    pub(crate) fn incr_count(&mut self, m: usize) {
        self.counts[m] += 1;
        if m >= 2 {
            self.counts[1] += 1
        }
    }

    #[inline]
    pub(crate) fn add_to_pileup(&mut self, ix: usize, e: PileupEntry) {
        self.pileup.as_mut().expect("Missing pileup")[ix] = e;
    }

    #[inline]
    pub(super) fn add_counts(&mut self, other: &Self) {
        for (p1, p2) in self.counts.iter_mut().zip(other.counts.iter()) {
            *p1 += *p2
        }
    }

    #[inline]
    pub(super) fn merge_pileups(&mut self, other: &Self) {
        if let Some(v1) = self.pileup.as_mut() {
            let v2 = other.pileup.as_deref().expect("Expected to find pileup");
            if v2.len() > v1.len() {
                v1.resize(v2.len(), PileupEntry::default())
            }
            for (p1, p2) in v1.iter_mut().zip(v2.iter()) {
                if p2.is_present() {
                    assert!(!p1.is_present(), "Clash when merging pileups");
                    *p1 = *p2
                }
            }
        } else {
            assert!(other.pileup.is_none(), "Didn't expect to find pileup");
        }
    }

    pub(super) fn reserve_pileup(&mut self, ix: usize) {
        let v = self.pileup.as_mut().expect("Missing pileup vector");

        let uncalled = PileupEntry::new_missing();

        if v.len() < ix + 1 {
            v.resize(ix + 1, PileupEntry::default())
        }
        v[ix] = uncalled;
    }

    #[inline]
    pub(crate) fn add_offset(&mut self, delta: u32) {
        self.offset += delta
    }

    #[inline]
    pub(crate) fn counts(&self) -> &[u32; 4] {
        &self.counts
    }

    #[inline]
    pub(crate) fn strand(&self) -> Strand {
        self.strand
    }

    #[inline]
    pub(crate) fn context(&self) -> CContext {
        self.context
    }

    #[inline]
    pub(crate) fn offset(&self) -> u32 {
        self.offset
    }

    #[inline]
    pub(crate) fn set_offset(&mut self, x: u32) {
        self.offset = x
    }

    #[inline]
    pub(crate) fn pileup(&self) -> Option<&[PileupEntry]> {
        self.pileup.as_deref()
    }
}
