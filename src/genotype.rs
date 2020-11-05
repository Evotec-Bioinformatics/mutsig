/// Resemble the genotype of a single sample. The genotype is a
/// list of allele-indizes that match the alleles given in the
/// VCF record.
///
/// Since we are only interested in the number of different alleles,
/// and not the haplotype, they are stored in a sorted vector.
pub struct Genotype {
    inner: Vec<Option<u8>>,
}

impl Genotype {
    /// Iterate the genotypes.
    pub fn iter(&self) -> GenotypeAlleleIterator {
        GenotypeAlleleIterator {
            inner: self.inner.iter().filter_map(|v| v.clone()).collect(),
            pos: 0,
        }
    }
}

use std::fmt;
impl fmt::Debug for Genotype {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.inner[0] {
            Some(i) => write!(f, "Genotype({}", i)?,
            None => write!(f, "Genotype(-")?,
        }
        for i in 1..self.inner.len() {
            match self.inner[i] {
                Some(g) => write!(f, "/{}", g)?,
                None => write!(f, "/-")?,
            }
        }
        write!(f, ")")
    }
}

impl From<rust_htslib::bcf::record::Genotype> for Genotype {
    fn from(gt: rust_htslib::bcf::record::Genotype) -> Genotype {
        // Reformat the rust-htslib allele indizes
        let mut inner: Vec<Option<u8>> = gt
            .iter()
            .map(|gta| match gta {
                rust_htslib::bcf::record::GenotypeAllele::Unphased(i) => Some(*i as u8),
                rust_htslib::bcf::record::GenotypeAllele::Phased(i) => Some(*i as u8),
                _ => None,
            })
            .collect();

        // Sort the vector (remember, we are only interested in the counts and not the haplotypes)
        inner.sort();
        Genotype { inner: inner }
    }
}

use std::cmp::*;
impl PartialEq<Genotype> for Genotype {
    fn eq(&self, other: &Genotype) -> bool {
        if self.inner.len() != other.inner.len() {
            return false;
        }

        for i in 0..self.inner.len() {
            if self.inner[i].is_none() || other.inner[i].is_none() {
                return true;
            }
            if self.inner[i] != other.inner[i] {
                return false;
            }
        }

        return true;
    }
}

/// Simple iterator of the genotypes of one sample
pub struct GenotypeAlleleIterator {
    inner: Vec<u8>,
    pos: usize,
}

impl Iterator for GenotypeAlleleIterator {
    type Item = u8;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.inner.len() {
            None
        } else {
            self.pos += 1;
            Some(self.inner[self.pos - 1])
        }
    }
}
