# MutSig

A slim tool for counting the mutational signature from a VCF file.

## Usage

Count number of single-nucleotide-variations in a single- or multi-sample VCF file:
```bash
mutsig my_sample.vcf.gz reference_genome.fa.gz > singlets.txt
``` 
The output will be a simple matrix in tab-separated text format containing the 
variant in the first column followed by the samples in the 

To count single-nucleotide variations in triplets, specify the window-size (i.e., number of bases up- and downstream to consider).
```bash
mutsig my_sample.vcf.gz reference_genome.fa.gz -w 1 > triplets.txt
``` 

If you have a multi-sample VCF you may want to ignore the position which are homogeneous in all
samples:
```bash
mutsig my_sample.vcf.gz reference_genome.fa.gz -i > non_homogeneous_singlets.txt
``` 

## Installation

### Cargo

If you have the Rust toolchain installed, it is as simple as

```bash
cargo install
``` 

### Conda

Planned

## Author

- Manuel Landesfeind <manuel.landesfeind@evotec.com>

