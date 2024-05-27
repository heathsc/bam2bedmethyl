use crate::read::pileup::PileupEntry;

pub(crate) struct CountBlock {
    idx: usize,
    tid: usize,
    start: usize,
    cpg_sites: Vec<CpG>,
}

impl CountBlock {
    /// Initialize a new count block.  Make a list of CpG sites in the reference sequence.   Note
    /// that rf has the reference only for the region covered by the block.
    pub(crate) fn new(idx: usize, tid: usize, start: usize, rf: &[u8], gen_pileup: bool) -> Self {
        let mut cpg_sites = Vec::new();

        // Unlikely to happen, but just in case...
        assert!(rf.len() < u32::MAX as usize, "Region too large");

        // Go through reference to find CpG sites
        for (i, c) in rf.windows(2).enumerate() {
            if c[0].to_ascii_uppercase() == b'C' && c[1].to_ascii_uppercase() == b'G' {
                cpg_sites.push(CpG::new(i as u32, gen_pileup))
            }
        }

        Self {
            idx,
            tid,
            start,
            cpg_sites,
        }
    }

    pub(crate) fn init_pileup_index(&mut self, ix: usize, x: usize, y: usize, reverse: bool) {
        assert!(x >= self.start && y >= x, "Illegal read coordinates");
        let x1 = x - self.start;
        let y1 = y - self.start;
        let i = self
            .cpg_sites
            .binary_search_by_key(&x1, |c| c.offset as usize)
            .unwrap_or_else(|i| i);
        let j = self
            .cpg_sites
            .binary_search_by_key(&y1, |c| c.offset as usize)
            .map(|i| i + 1)
            .unwrap_or_else(|i| i);
        for c in self.cpg_sites[i..j].iter_mut() {
            c.reserve_pileup(ix, reverse)
        }
    }

    #[inline]
    pub(crate) fn start(&self) -> usize {
        self.start
    }

    #[inline]
    pub(crate) fn idx(&self) -> usize {
        self.idx
    }

    #[inline]
    pub(crate) fn tid(&self) -> usize {
        self.tid
    }

    #[inline]
    pub(crate) fn cpg_sites(&self) -> &[CpG] {
        &self.cpg_sites
    }

    #[inline]
    pub(crate) fn cpg_sites_mut(&mut self) -> &mut Vec<CpG> {
        &mut self.cpg_sites
    }

    pub(crate) fn find_site(&mut self, x: u32) -> Option<&mut CpG> {
        self.cpg_sites
            .binary_search_by_key(&x, |c| c.offset)
            .map(|i| &mut self.cpg_sites[i])
            .ok()
    }

    #[inline]
    pub(crate) fn end(&self) -> usize {
        self.start + self.cpg_sites.last().map_or(0, |c| c.offset as usize)
    }

    pub(crate) fn add_counts(&mut self, other: &mut Self) {
        assert!(self.start >= other.start);
        // Only add counts if overlapping
        if self.start <= other.end() && !self.cpg_sites.is_empty() {
            let delta = (self.start - other.start) as u32;
            let first_x = self.cpg_sites[0].offset + delta;
            let ix = other
                .cpg_sites
                .binary_search_by_key(&first_x, |c| c.offset)
                .expect("CpG site not found!");
            for (c1, c2) in other.cpg_sites[ix..].iter().zip(self.cpg_sites.iter_mut()) {
                c2.add_counts(c1);
                c2.merge_pileups(c1);
            }
            if other.end() > self.end() {
                let last_x = self.cpg_sites.last().map_or(0, |c| c.offset) + delta;
                let ix = other
                    .cpg_sites
                    .binary_search_by_key(&last_x, |c| c.offset)
                    .expect("CpG site not found!");
                for mut c in other.cpg_sites.drain(ix + 1..) {
                    c.offset = c.offset.checked_sub(delta).expect("Illegal offset");
                    self.cpg_sites.push(c)
                }
            }
        }
    }
}
/// The counts are:
///   0 - not modified
///   1 - 5mC or 5hmC
///   2 - 5mC
///   3 - 5hmC
pub(super) struct CpG {
    offset: u32,
    fwd_counts: [u32; 4],
    rev_counts: [u32; 4],
    pileup: Option<Vec<PileupEntry>>,
}

impl CpG {
    #[inline]
    fn new(offset: u32, gen_pileup: bool) -> Self {
        Self {
            offset,
            fwd_counts: [0; 4],
            rev_counts: [0; 4],
            pileup: if gen_pileup { Some(Vec::new()) } else { None },
        }
    }

    #[inline]
    pub(super) fn incr_count(&mut self, m: usize, reverse: bool) {
        let ct = if reverse {
            &mut self.rev_counts
        } else {
            &mut self.fwd_counts
        };
        ct[m] += 1;
        if m >= 2 {
            ct[1] += 1
        }
    }

    #[inline]
    pub(super) fn add_to_pileup(&mut self, ix: usize, m: PileupEntry) {
        self.pileup.as_mut().expect("Missing pileup")[ix] = m;
    }

    #[inline]
    fn add_counts(&mut self, other: &Self) {
        for (p1, p2) in self.fwd_counts.iter_mut().zip(other.fwd_counts.iter()) {
            *p1 += *p2
        }
        for (p1, p2) in self.rev_counts.iter_mut().zip(other.rev_counts.iter()) {
            *p1 += *p2
        }
    }

    #[inline]
    fn merge_pileups(&mut self, other: &Self) {
        if let Some(v1) = self.pileup.as_mut() {
            let v2 = other.pileup.as_deref().expect("Expected to find pileup");
            if v2.len() > v1.len() {
                v1.resize(v2.len(), PileupEntry::default())
            }
            for (p1, p2) in v1.iter_mut().zip(v2.iter()) {
                match (*p1, *p2) {
                    (_, PileupEntry::NotPresent) => {}
                    (PileupEntry::NotPresent, e) => *p1 = e,
                    _ => panic!("Clash when merging pileups"),
                }
            }
        } else {
            assert!(other.pileup.is_none(), "Didn't expect to find pileup");
        }
    }

    fn reserve_pileup(&mut self, ix: usize, reverse: bool) {
        let v = self.pileup.as_mut().expect("Missing pileup vector");

        let uncalled = if reverse {
            PileupEntry::UncalledRev
        } else {
            PileupEntry::UncalledFwd
        };

        if v.len() < ix + 1 {
            v.resize(ix + 1, PileupEntry::default())
        }
        v[ix] = uncalled;
    }

    #[inline]
    pub(super) fn add_offset(&mut self, delta: u32) {
        self.offset += delta
    }

    #[inline]
    pub(super) fn fwd_counts(&self) -> &[u32; 4] {
        &self.fwd_counts
    }

    #[inline]
    pub(super) fn rev_counts(&self) -> &[u32; 4] {
        &self.rev_counts
    }

    #[inline]
    pub(super) fn offset(&self) -> u32 {
        self.offset
    }

    #[inline]
    pub(super) fn pileup(&self) -> Option<&[PileupEntry]> {
        self.pileup.as_deref()
    }
}
