/// The result matrix that contains the results
/// by means of counts of variants per sample.
pub struct ResultMatrix {
    n_samples: usize,
    inner: Vec<u32>,
}

impl ResultMatrix {
    /// Create a new matrix containing data for `n_variants` and `n_samples`.
    pub fn new(n_variants: usize, n_samples: usize) -> Self {
        let n_total = n_variants * n_samples;
        let v = (0..n_total).map(|_| 0).collect();
        ResultMatrix {
            n_samples: n_samples,
            inner: v,
        }
    }

    /// Calculate the index position in `inner` for given variant index `vidx and sample index `sidx`.
    fn index(&self, vidx: usize, sidx: usize) -> usize {
        vidx * self.n_samples + sidx
    }

    /// Increment the count for variant at `vidx` and samples at `sidx` by one.
    pub fn increment(&mut self, vidx: usize, sidx: usize) {
        let idx = self.index(vidx, sidx);
        self.inner[idx] += 1;
    }

    /// Return the count for variant at `vidx` and sample at `sidx`.
    pub fn get(&mut self, vidx: usize, sidx: usize) -> u32 {
        self.inner[self.index(vidx, sidx)]
    }
}
