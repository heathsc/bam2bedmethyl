#[derive(Copy, Clone, Debug)]
pub struct OptIndex {
    inner: i64,
}

impl Default for OptIndex {
    fn default() -> Self {
        Self { inner: -1 }
    }
}

impl OptIndex {
    #[inline]
    pub fn get(&self) -> Option<usize> {
        if self.inner < 0 {
            None
        } else {
            Some(self.inner as usize)
        }
    }
    #[inline]
    pub fn set(&mut self, x: Option<usize>) {
        if let Some(y) = x {
            self.inner = y.try_into().expect("Illegal index value")
        } else {
            self.inner = -1
        }
    }
}
