use super::{
    Config,
    cytosine::{CContext, Cytosine, Strand},
};

pub(crate) struct CountBlock {
    idx: usize,
    tid: usize,
    start: usize,
    c_sites: Vec<Cytosine>,
}

impl CountBlock {
    /// Initialize a new count block.  Make a list of CpG sites in the reference sequence.   Note
    /// that rf has the reference only for the region covered by the block.
    pub(crate) fn new(
        idx: usize,
        tid: usize,
        start: usize,
        rf: &[u8],
        cfg: &Config,
        trim: usize,
    ) -> Self {
        let mut c_sites = Vec::new();
        let gen_pileup = cfg.pileup();
        let non_cpg = cfg.non_cpg();

        let l = rf.len();
        // Unlikely to happen, but just in case...
        assert!(l < u32::MAX as usize, "Region too large");

        let mut add_site = |i: usize, s: Strand, c: CContext, p: bool| {
            if i > trim && i < l - trim {
                c_sites.push(Cytosine::new(i as u32, s, c, p))
            }
        };

        // Go through reference to find C sites
        if non_cpg {
            for (i, c) in rf.windows(3).enumerate() {
                if c[0].eq_ignore_ascii_case(&b'C') {
                    if c[1].eq_ignore_ascii_case(&b'G') {
                        add_site(i, Strand::Forward, CContext::Cg, gen_pileup);
                        add_site(i + 1, Strand::Reverse, CContext::Cg, gen_pileup);
                        if c[2].eq_ignore_ascii_case(&b'G') {
                            add_site(i + 2, Strand::Reverse, CContext::Chg, false)
                        }
                    } else if c[2].eq_ignore_ascii_case(&b'G') {
                        add_site(i, Strand::Forward, CContext::Chg, false);
                        if !c[1].eq_ignore_ascii_case(&b'C') {
                            add_site(i + 2, Strand::Reverse, CContext::Chg, false)
                        }
                    } else {
                        add_site(i, Strand::Forward, CContext::Chh, false)
                    }
                } else if c[2].eq_ignore_ascii_case(&b'G')
                    && !c[1].eq_ignore_ascii_case(&b'C')
                    && !c[0].eq_ignore_ascii_case(&b'C')
                {
                    add_site(i + 2, Strand::Reverse, CContext::Chh, false)
                }
            }
        } else {
            // CpGs only
            for (i, c) in rf.windows(2).enumerate() {
                if c[0].eq_ignore_ascii_case(&b'C') && c[1].eq_ignore_ascii_case(&b'G') {
                    add_site(i, Strand::Forward, CContext::Cg, gen_pileup);
                    add_site(i + 1, Strand::Reverse, CContext::Cg, gen_pileup);
                }
            }
        }

        if non_cpg {
            c_sites.sort_unstable_by_key(|c| c.offset());
            for c in c_sites.windows(2) {
                if c[0].offset() == c[1].offset() {
                    eprintln!("OOOOOK! {c:?}")
                }
            }
        }

        Self {
            idx,
            tid,
            start,
            c_sites,
        }
    }

    pub(crate) fn init_pileup_index(&mut self, ix: usize, x: usize, y: usize) {
        assert!(x >= self.start && y >= x, "Illegal read coordinates");
        let x1 = x - self.start;
        let y1 = y - self.start;
        let mut i = self
            .c_sites
            .binary_search_by_key(&x1, |c| c.offset() as usize)
            .unwrap_or_else(|i| i);

        while self.c_sites[i].context() != CContext::Cg {
            if i > 0 {
                i -= 1;
            } else {
                break;
            }
        }
        if i > 0 && self.c_sites[i].strand() == Strand::Reverse {
            i -= 1;
        }

        let mut j = self
            .c_sites
            .binary_search_by_key(&y1, |c| c.offset() as usize)
            .unwrap_or_else(|i| i);

        let l = self.c_sites.len();
        while self.c_sites[j].context() != CContext::Cg {
            if j + 1 < l { j += 1 } else { break }
        }
        for c in self.c_sites[i..=j].iter_mut() {
            c.reserve_pileup(ix)
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
    pub(crate) fn c_sites(&self) -> &[Cytosine] {
        &self.c_sites
    }

    #[inline]
    pub(crate) fn c_sites_mut(&mut self) -> &mut Vec<Cytosine> {
        &mut self.c_sites
    }

    pub(crate) fn find_site(&mut self, x: u32, r: bool) -> Option<&mut Cytosine> {
        let s = if r { Strand::Reverse } else { Strand::Forward };

        self.c_sites
            .binary_search_by_key(&x, |c| c.offset())
            .map(|i| &mut self.c_sites[i])
            .ok()
            .and_then(|c| if c.strand() == s { Some(c) } else { None })
    }

    #[inline]
    pub(crate) fn end(&self) -> usize {
        self.start + self.c_sites.last().map_or(0, |c| c.offset() as usize)
    }

    pub(crate) fn add_counts(&mut self, other: &mut Self) {
        assert!(self.start >= other.start);
        // Only add counts if overlapping
        if self.start <= other.end() && !self.c_sites.is_empty() {
            let delta = (self.start - other.start) as u32;
            let first_x = self.c_sites[0].offset() + delta;
            let ix = other
                .c_sites
                .binary_search_by_key(&first_x, |c| c.offset())
                .expect("Cytosine site not found!");
            for (c1, c2) in other.c_sites[ix..].iter().zip(self.c_sites.iter_mut()) {
                c2.add_counts(c1);
                c2.merge_pileups(c1);
            }
            if other.end() > self.end() {
                let last_x = self.c_sites.last().map_or(0, |c| c.offset()) + delta;
                let ix = other
                    .c_sites
                    .binary_search_by_key(&last_x, |c| c.offset())
                    .expect("Cytosine site not found!");
                for mut c in other.c_sites.drain(ix + 1..) {
                    c.set_offset(c.offset().checked_sub(delta).expect("Illegal offset"));
                    self.c_sites.push(c)
                }
            }
        }
    }
}
