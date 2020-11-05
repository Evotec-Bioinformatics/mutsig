use rust_htslib::faidx;
use std::path::Path;

/// Struct for fetching data from a faidx'ed FASTQ file. Automatically retrieves
/// also the surrounding bases as requested per `window_size`.
pub struct Reference {
    inner: faidx::Reader,
    window: u8,
}

impl Reference {
    /// Create a new reference setup using the FASTA file located at
    /// `path` and a window size of `window_size` bases.
    ///
    /// # Arguments
    /// * `window_size` the number of bases to retrieve up- and down-stream of the requested position. A `window_size` of 1 will return triplets in `fetch()`.
    pub fn new<P: AsRef<Path>>(path: P, window_size: u8) -> Result<Self, String> {
        let inner = match faidx::Reader::from_path(path.as_ref()) {
            Err(e) => {
                return Err(format!(
                    "Can not open '{}': {}",
                    path.as_ref().to_str().unwrap(),
                    e
                ))
            }
            Ok(i) => i,
        };

        Ok(Reference {
            inner: inner,
            window: window_size,
        })
    }

    /// Get the reference sequence at a given position (0-based offset). If a window-size was given during
    /// creation of the reference, then that number of bases before and after position are extracted too.
    pub fn fetch<N: AsRef<str>>(&self, name: N, position: i64) -> Result<String, String> {
        if self.window as i64 > position {
            return Err(format!(
                "Can not fetch window {} before {}",
                self.window, position
            ));
        }
        let start = position - self.window as i64;
        let end = position + self.window as i64;

        Ok(self
            .inner
            .fetch_seq_string(name, start as usize, end as usize)
            .unwrap()
            .to_uppercase()
            .to_string())
    }

    /// Retrieve the window size
    pub fn window_size(&self) -> u8 {
        self.window
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build(window_size: u8) -> Reference {
        Reference::new(
            format!("{}/testdata/ex2.fa", env!("CARGO_MANIFEST_DIR")),
            window_size,
        )
        .ok()
        .unwrap()
    }

    #[test]
    fn test_singlet_start_of_chromosome() {
        let r = build(0).fetch("1", 0);
        assert_eq!(r, Ok("T".to_owned()))
    }

    #[test]
    fn test_singlet_end_of_chromosome() {
        let r = build(0).fetch("1", 5);
        assert_eq!(r, Ok("A".to_owned()))
    }

    #[test]
    fn test_triplet_end_of_chromosome() {
        let r = build(1).fetch("1", 5);
        assert_eq!(r, Ok("GA".to_owned()))
    }
}
