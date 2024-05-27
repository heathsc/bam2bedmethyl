use rs_htslib::sam::{ModIter, SeqBase};

// Combine mod probabilities for the same position (i.e., for 5hmC and 5mC)
pub(super) struct MethItr<'a, 'b> {
    mod_iter: ModIter<'a, 'b>,
    seq_ix: usize,
    finished: bool,
}

impl<'a, 'b> MethItr<'a, 'b> {
    pub(super) fn new(mod_iter: ModIter<'a, 'b>) -> Self {
        Self {
            mod_iter,
            seq_ix: 0,
            finished: false,
        }
    }
}

impl<'a, 'b> Iterator for MethItr<'a, 'b> {
    type Item = (usize, SeqBase, Option<(u8, u8, bool)>);

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            None
        } else if let Some(d) = self.mod_iter.next_pos() {
            let i = self.seq_ix;
            let b = d.seq_base();
            let v = d.data();
            let p = if v.is_empty() {
                None
            } else {
                let (tp_m, tp_h) =
                    v.iter()
                        .fold((0, 0), |(t_m, t_h), (p, m)| match m.base_mod_code() {
                            Some(b'm') => (t_m + (*p as usize), t_h),
                            Some(b'h') => (t_m, t_h + (*p as usize)),
                            _ => (t_m, t_h),
                        });
                Some((
                    tp_m.min(255) as u8,
                    tp_h.min(255) as u8,
                    v[0].1.reverse_strand(),
                ))
            };
            self.seq_ix += 1;
            Some((i, b, p))
        } else {
            self.finished = true;
            None
        }
    }
}
