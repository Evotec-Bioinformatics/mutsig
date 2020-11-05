extern crate pretty_env_logger;
#[macro_use]
extern crate log;
extern crate rust_htslib;
use rust_htslib::bcf::Read;
use std::collections::BTreeMap;
mod genotype;
mod reference;
mod result;
mod signature;

fn main() -> Result<(), String> {
    #[cfg(debug_assertions)]
    {
        if std::env::var("RUST_LOG").is_err() {
            std::env::set_var("RUST_LOG", "trace");
        }
    }

    pretty_env_logger::init();
    let matches = clap::App::new(env!("CARGO_PKG_NAME"))
        .version(env!("CARGO_PKG_VERSION"))
        .author(env!("CARGO_PKG_AUTHORS"))
        .arg(
            clap::Arg::with_name("VCF")
                .help("Sets the input VCF file to use")
                .required(true)
                .index(1),
        )
        .arg(
            clap::Arg::with_name("REFERENCE")
                .help("Sets the input reference FASTQ file (must be indexed with faidx)")
                .required(true)
                .index(2),
        )
        .arg(
            clap::Arg::with_name("samples")
                .short("s")
                .long("samples")
                .help("Include this sample in the analyzes (defaults to all), can be specified multiple times")
								.takes_value(true)
								.value_name("SAMPLE")
								.multiple(true)
        )
        .arg(
            clap::Arg::with_name("ignore-homogeneous")
                .short("i")
                .long("ignore-homogeneous")
                .help("Ignore sites where all samples have the same allele count")
        )
        .arg(
            clap::Arg::with_name("window")
                .short("w")
                .long("bases-window")
                .help("The number of bases to consider up and downstream of the mutation position")
                .value_name("BASES")
                .takes_value(true),
        )
        .get_matches();
    info!(
        "Started {} v{}",
        env!("CARGO_PKG_NAME"),
        env!("CARGO_PKG_VERSION")
    );

    let ignore_homogeneous_sites = matches.occurrences_of("ignore-homogeneous") > 0;

    // Window size parameter
    let window_size = match matches.value_of("window") {
        Some(v) => match v.parse::<u8>() {
            Err(e) => return Err(format!("Invalid window-parameter '{}': {}", v, e)),
            Ok(w) => w,
        },
        None => 0,
    };

    // Load access to the reference
    let reference = match matches.value_of("REFERENCE") {
        None => return Err("Require 'REFERENCE' file name".into()),
        Some(path) => {
            info!(
                "Using reference from {} with window size of {}",
                path, window_size
            );
            match reference::Reference::new(path, window_size) {
                Ok(r) => r,
                Err(e) => return Err(format!("Can not open reference '{}': {}", path, e)),
            }
        }
    };

    // Open the VCF file
    let mut variants = match matches.value_of("VCF") {
        None => return Err("Require 'VCF' file name".into()),
        Some(p) => match rust_htslib::bcf::Reader::from_path(p) {
            Err(e) => return Err(format!("Can not open VCF file '{}': {}", p, e)),
            Ok(v) => v,
        },
    };

    // Fetch information about the contigs
    let header = variants.header();
    let contigs: BTreeMap<u32, String> = (0..header.contig_count())
        .map(|rid| {
            (
                rid,
                std::str::from_utf8(header.rid2name(rid).ok().unwrap())
                    .unwrap()
                    .to_owned(),
            )
        })
        .collect();

    // Fetch information about the samples.
    let sample_names: Vec<String> = header
        .samples()
        .iter()
        .map(|u| std::str::from_utf8(u).unwrap().to_owned())
        .collect();
    // In case a list of samples to analyze is requested, find their indizes
    // and memorize them.
    let bcf_sample_indizes: Vec<usize> = match matches.values_of("samples") {
        None => (0..sample_names.len()).collect(),
        Some(values) => {
            let mut s = Vec::new();
            for v in values {
                s.push(index_of(&sample_names, v)?)
            }
            s
        }
    };

    let n_samples = bcf_sample_indizes.len();
    debug!(
        "Processing {} samples ({:?}) from: {:?}",
        n_samples, bcf_sample_indizes, sample_names
    );

    // We can only ignore homgeneous sites if we have more than one sample
    if ignore_homogeneous_sites && n_samples < 2 {
        return Err(
            "Found only one sample but were told to ignor homgeneous sites - this is not possible"
                .to_owned(),
        );
    }

    // Build a list of all signatures
    let signatures = signature::Signatures::new(window_size.into());
    let n_variants = signatures.len();
    debug!("Found a total of {} signature variants", n_variants);

    // Initialize the result matrix
    let mut results = result::ResultMatrix::new(n_variants, n_samples);

    // Iterate the codonds
    for res_record in variants.records() {
        let mut record = match res_record {
            Ok(r) => r,
            Err(e) => return Err(format!("Can not retrieve next VCF record: {}", e)),
        };

        // Fetch all the alleles
        let alleles = match alternative_alleles_from_record(&record, &contigs, &reference) {
            AlleleRecordStatus::Ok(a) => a,
            AlleleRecordStatus::Err(e) => return Err(e),
            AlleleRecordStatus::Ignore(e) => {
                trace!("{}", e);
                continue;
            }
            AlleleRecordStatus::Issue(e) => {
                warn!("{}", e);
                continue;
            }
        };
        debug!("Found alleles: {:?}", alleles);

        // Match the allele(-indize)s into the signature_indizes
        let signature_indizes: Vec<usize> = alleles
            .iter()
            .map(|a| signatures.index_of(a).unwrap())
            .collect();
        debug!("Found signature indizes: {:?}", signature_indizes);

        // Extract the genotypes from the record in the order of our
        // expected/wanted samples and re-encode them as our genotype struct
        let bcf_gts = record.genotypes().unwrap();
        let gts = bcf_sample_indizes
            .iter()
            .map(|sample_index| genotype::Genotype::from(bcf_gts.get(*sample_index)))
            .collect();
        trace!("Found genotypes: {:?}", gts);

        // If all sites should be counted or there is variance in the genotypes
        if !ignore_homogeneous_sites || is_varying_position(&gts) {
            // for each sample
            for sample_index in 0..n_samples {
                // for each allele of that sample
                for allele_index in gts[sample_index].iter() {
                    // if it is not the reference
                    if allele_index > 0 {
                        // get the signature and increment it
                        let sig_index = signature_indizes[allele_index as usize - 1];
                        results.increment(sig_index, sample_index)
                    }
                }
            }
        }
    }

    // Identify the signatures that we want to report
    let forwards: Vec<signature::Signature> = signatures
        .signatures()
        .iter()
        .filter(|s| s.is_forward_signature())
        .map(|c| c.clone())
        .collect();

    // Print header
    print!("Variant");
    for sidx in bcf_sample_indizes {
        print!("\t{}", sample_names[sidx]);
    }
    println!("");

    // Print the results
    for v in 0..forwards.len() {
        let signature = &forwards[v];
        let signature_index = signatures.index_of(signature).unwrap();
        print!("{}", signature);
        for s in 0..n_samples {
            print!("\t{}", results.get(signature_index, s));
        }
        println!("");
    }

    Ok(())
}

enum AlleleRecordStatus {
    Err(String),
    Issue(String),
    Ignore(String),
    Ok(Vec<signature::Signature>),
}

/// Extract the alternative alleles from a VCF record.
/// If succesful, a list of signatures resembeling the codon-allele combinations is returned.
fn alternative_alleles_from_record(
    record: &rust_htslib::bcf::Record,
    contigs: &BTreeMap<u32, String>,
    reference: &reference::Reference,
) -> AlleleRecordStatus {
    let mut alleles = Vec::new();

    // Identify contig as string
    let contig = match contigs.get(&record.rid().unwrap()) {
        Some(c) => c,
        None => {
            return AlleleRecordStatus::Err(format!(
                "Can not find contig name for template-id {}",
                record.rid().unwrap()
            ))
        }
    };
    let position = record.pos() as usize;

    // Iterator on all alleles
    let record_alleles = record.alleles();
    let mut allele_iter = record_alleles.iter();
    // Expect the first allele to be the reference allele
    let reference_allele = std::str::from_utf8(allele_iter.next().unwrap())
        .unwrap()
        .to_uppercase();

    // Ignore deletion events
    if reference_allele.len() > 1 {
        return AlleleRecordStatus::Ignore(format!(
            "Ignoring non-SNV variant at position {}:{}",
            contig,
            position + 1
        ));
    }
    let reference_nucleotide = reference_allele.as_bytes()[0] as char;
    let codon = match reference.fetch(contig, position as i64) {
        Ok(s) => s,
        Err(e) => {
            return AlleleRecordStatus::Err(format!(
                "Can not fetch codon at position {}:{}: {}",
                contig, position, e
            ))
        }
    };
    // Check that the codon is ACGT only
    if codon
        .chars()
        .any(|c| c != 'A' && c != 'C' && c != 'G' && c != 'T')
    {
        return AlleleRecordStatus::Ignore(format!(
            "Ignoring codon with non-standard nucleotide at position {}:{}: {}",
            contig,
            position + 1,
            codon
        ));
    }

    // Check that the nucleotide at the reference-position in the codon matches the reference-allele from the record
    if codon.as_bytes()[reference.window_size() as usize] as char != reference_nucleotide {
        return AlleleRecordStatus::Issue(format!(
            "Loaded codon '{}' does not match to expected reference allele {}",
            codon, reference_nucleotide
        ));
    }
    // Ensure that no allele is a insertion
    for a in allele_iter {
        if a.len() != 1 {
            return AlleleRecordStatus::Ignore(format!(
                "Ignoring non-SNV variant at position {}:{}",
                contig,
                position + 1
            ));
        } else {
            // for SNPs, push the signature to the result list
            alleles.push(signature::Signature::new(
                &codon,
                reference_nucleotide,
                (a[0] as char).to_uppercase().next().unwrap(),
            ))
        }
    }

    if alleles.len() < 1 {
        return AlleleRecordStatus::Ignore(format!(
            "Ignoring no-alternative variant at position {}:{}",
            contig,
            position + 1
        ));
    }

    AlleleRecordStatus::Ok(alleles)
}

/// Helper function to check if there is variation in the genotypes
fn is_varying_position(gts: &Vec<genotype::Genotype>) -> bool {
    for i in 1..gts.len() {
        if gts[0] != gts[i] {
            return true;
        }
    }
    return false;
}

/// Helper function to find the index of a sample in the sample list
fn index_of(hay: &Vec<String>, needle: &str) -> Result<usize, String> {
    for i in 0..hay.len() {
        if needle == hay[i] {
            return Ok(i);
        }
    }

    let mut msg = format!("Can not find sample '{}' in list of: '{}'", needle, hay[0]);
    for h in &hay[1..] {
        msg.push_str(", '");
        msg.push_str(h);
        msg.push_str("'");
    }
    Err(msg)
}

#[cfg(test)]
mod tests {}
