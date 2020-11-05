use std::cmp;
use std::collections::BTreeMap;
use std::fmt;

#[derive(Clone)]
pub struct Signature {
    codon: String,
    reference: char,
    alternative: char,
}
impl fmt::Debug for Signature {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}>{}", self.codon, self.reference, self.alternative)
    }
}

impl Signature {
    pub fn new<S: AsRef<str>>(codon: S, reference: char, alternative: char) -> Signature {
        Signature {
            codon: codon.as_ref().to_owned(),
            reference: reference,
            alternative: alternative,
        }
    }

    pub fn is_forward_signature(&self) -> bool {
        self.reference == 'C' || self.reference == 'T'
    }
}

impl fmt::Display for Signature {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}>{}", self.codon, self.alternative)
    }
}

impl cmp::Ord for Signature {
    fn cmp(&self, other: &Self) -> cmp::Ordering {
        let c = self.codon.cmp(&other.codon);
        if c != cmp::Ordering::Equal {
            c
        } else {
            self.alternative.cmp(&other.alternative)
        }
    }
}

impl cmp::PartialOrd for Signature {
    fn partial_cmp(&self, other: &Self) -> Option<cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl cmp::PartialEq for Signature {
    fn eq(&self, other: &Self) -> bool {
        self.codon == other.codon && self.alternative == other.alternative
    }
}
impl cmp::Eq for Signature {}

#[derive(Clone)]
pub struct Signatures {
    db: BTreeMap<Signature, usize>,
}

impl Signatures {
    pub fn new(window: usize) -> Signatures {
        let s = Signatures {
            db: build_signatures(window),
        };
        trace!("Build signature database: {:?}", s);
        s
    }

    /// Returns the index of the signature
    pub fn index_of(&self, sig: &Signature) -> Option<usize> {
        match self.db.get(sig) {
            Some(i) => Some(i.clone()),
            None => None,
        }
    }

    /// Retrieve the number of signatures
    pub fn len(&self) -> usize {
        self.db.len()
    }

    pub fn signatures(&self) -> Vec<Signature> {
        self.db.keys().map(|i| i.clone()).collect()
    }
}

impl fmt::Debug for Signatures {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Signatures {{")?;
        for (k, v) in &self.db {
            write!(f, "{:?} = {}, ", k, v)?;
        }
        write!(f, "}}")
    }
}

fn build_codons(is_forward: bool, mut prior: usize, mut after: usize) -> Vec<String> {
    let mut codons: Vec<String> = vec!["".to_owned()];

    while prior >= 1 {
        let mut sig2: Vec<String> = Vec::new();
        for c in vec!['A', 'C', 'G', 'T'] {
            for s in &codons {
                let mut s2 = s.clone();
                s2.push(c);
                sig2.push(s2);
            }
        }
        codons = sig2;
        prior -= 1;
    }

    let reference_nucleotides = if is_forward {
        vec!['C', 'T']
    } else {
        vec!['A', 'G']
    };

    let mut sig2: Vec<String> = Vec::new();
    for c in reference_nucleotides {
        for s in &codons {
            let mut s2 = s.clone();
            s2.push(c);
            sig2.push(s2);
        }
    }
    codons = sig2;

    while after >= 1 {
        let mut sig2: Vec<String> = Vec::new();
        for c in vec!['A', 'C', 'G', 'T'] {
            for s in &codons {
                let mut s2 = s.clone();
                s2.push(c);
                sig2.push(s2);
            }
        }
        codons = sig2;
        after -= 1;
    }

    codons
}

pub fn build_signatures(window: usize) -> BTreeMap<Signature, usize> {
    let mut signatures = BTreeMap::new();

    let mut idx = 0;
    for cod in build_codons(true, window, window) {
        let cn = cod.as_bytes()[window] as char;
        for n in vec!['A', 'C', 'G', 'T'] {
            if n != cn {
                signatures.insert(
                    Signature {
                        codon: cod.clone(),
                        reference: cn,
                        alternative: n,
                    },
                    idx,
                );
                idx += 1;
            }
        }
    }
    trace!("Built forward codons: {:?}", signatures);

    for cod in build_codons(false, window, window) {
        let cn = cod[window..(window + 1)].as_bytes()[0] as char;
        let cod2 = rev_comp(cod.chars());
        for n in vec!['A', 'C', 'G', 'T'] {
            if n != cn {
                let rev_sig = Signature::new(&cod2, cn, rev_comp_c(n));
                let idx = signatures.get(&rev_sig).unwrap().clone();
                signatures.insert(
                    Signature {
                        codon: cod.clone(),
                        reference: cn,
                        alternative: n,
                    },
                    idx,
                );
            }
        }
    }

    signatures
}

fn rev_comp<I: DoubleEndedIterator<Item = char>>(chars: I) -> String {
    chars.rev().map(rev_comp_c).collect::<String>()
}

fn rev_comp_c(n: char) -> char {
    match n {
        'A' => 'T',
        'C' => 'G',
        'G' => 'C',
        'T' => 'A',
        _ => panic!("Can not complement '{}'", n),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_rev_comp1() {
        assert_eq!(rev_comp("TGA".chars()), "TCA");
        assert_eq!(rev_comp("AGA".chars()), "TCT");
    }
}
