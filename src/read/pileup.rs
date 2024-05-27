use std::fmt::Formatter;
use std::{
    cmp::{Eq, Ord, Ordering, PartialEq, PartialOrd},
    collections::{BTreeSet, BinaryHeap},
    fmt,
};

struct PileupLease {
    index: usize,
    end_pos: usize,
}

impl Ord for PileupLease {
    fn cmp(&self, other: &Self) -> Ordering {
        other.end_pos.cmp(&self.end_pos)
    }
}

impl PartialOrd for PileupLease {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for PileupLease {
    fn eq(&self, other: &Self) -> bool {
        self.end_pos == other.end_pos
    }
}

impl Eq for PileupLease {}

#[derive(Default)]
pub(crate) struct Pileup {
    leases: BinaryHeap<PileupLease>,
    avail: BTreeSet<usize>,
}

impl Pileup {
    pub(crate) fn new() -> Self {
        Self::default()
    }

    pub(crate) fn new_index(&mut self, end_pos: usize) -> usize {
        let index = self.avail.pop_first().unwrap_or(self.leases.len());
        let pl = PileupLease { index, end_pos };
        self.leases.push(pl);
        index
    }

    pub(crate) fn trim(&mut self, start_pos: usize) {
        while let Some(pl) = self.leases.peek() {
            if pl.end_pos >= start_pos {
                break;
            }
            let ix = pl.index;
            self.avail.insert(ix);
            self.leases.pop();
        }
    }

    pub(crate) fn clear(&mut self) {
        self.leases.clear();
        self.avail.clear();
    }
}

#[derive(Debug, Default, Copy, Clone, Eq, PartialEq)]
pub(crate) enum PileupEntry {
    #[default]
    NotPresent = 0,
    CytosineFwd,
    MethCytosineFwd,
    HydroxyMethCytosineFwd,
    TotalMethFwd,
    UncalledFwd,
    CytosineRev,
    MethCytosineRev,
    HydroxyMethCytosineRev,
    TotalMethRev,
    UncalledRev,
}

const PE_OUTPUT: [char; 11] = [' ', 'C', 'M', 'H', 'T', '.', 'c', 'm', 'h', 't', ','];

impl fmt::Display for PileupEntry {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{}", PE_OUTPUT[*self as usize])
    }
}
