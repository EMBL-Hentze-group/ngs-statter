use crate::gff_parser::{self, gff3_parser};
use rust_htslib::bam::{self, IndexedReader, Read};
use std::collections::{HashMap, HashSet};

/// Computes the read distribution based on gene types from a BAM file and a gene map.
///
/// Given a BAM file, a gene map, and a minimum quality threshold, this function calculates
/// the distribution of reads for each gene type. It returns a HashMap where keys are gene types,
/// and values are inner HashMaps containing read length (u16) and corresponding read counts (u32).
///
/// # Arguments
///
/// * `bam_path` - The path to the BAM file.
/// * `gene_map` - A HashMap containing gene information, where keys are chromosome names and
///               values are vectors of `Gene` structs.
/// * `min_q` - The minimum quality threshold for considering a read.
///
/// # Returns
///
/// A HashMap containing gene types and their associated read distributions.
///
/// # Panics
///
/// Panics if there are issues reading the BAM file, parsing records, or if there are no common
/// chromosomes between the BAM file and the gene map.
///
/// # Examples
///
/// ```
/// use std::collections::{HashMap, HashSet};
/// use rust_htslib::bam::IndexedReader;
/// use your_crate_name::{compute_gene_type_read_dist, gff_parser::Gene};
///
/// let bam_path = "path/to/your.bam";
/// let gene_map: HashMap<String, Vec<Gene>> = /* initialize your gene map */;
/// let min_q = 20;
///
/// let result = compute_gene_type_read_dist(bam_path, gene_map, min_q);
/// ```
fn compute_gene_type_read_dist(
    bam_path: &str,
    gene_map: HashMap<String, Vec<gff_parser::Gene>>,
    min_q: u8,
    ignore_duplicate: bool,
) -> HashMap<String, HashMap<u16, u32>> {
    let mut bam: IndexedReader = match IndexedReader::from_path(bam_path) {
        Ok(abam) => abam,
        Err(e) => panic!("Cannot read {}: {}", bam_path, e.to_string()),
    };
    // get chromosome info from the bam header
    let header = bam::Header::from_template(bam.header()).to_hashmap();
    let chrom_data = match header.get("SQ") {
        Some(chrom) => chrom,
        None => panic!(
            "{:?} bam file does not contain chromosome info in header. Cannot parse this file!",
            bam_path
        ),
    };
    // common chromosomes in this bam file and input gff3 file
    let mut common_chroms: HashSet<String> = HashSet::new();
    for cd in chrom_data {
        let chrom = cd.get("SN").unwrap().to_owned();
        if gene_map.contains_key(&chrom) {
            common_chroms.insert(chrom);
        }
    }
    if common_chroms.len() == 0 {
        panic!("There are no common chromosomes in the given bam file and gff3 file. Check your inputs!");
    }
    let mut gene_type_read_dist: HashMap<String, HashMap<u16, u32>> = HashMap::new();
    for chrom in common_chroms {
        for gene in gene_map.get(&chrom).unwrap() {
            let gene_strand: bool = if gene.strand == "-" { true } else { false };
            if let Ok(_bfetch) = bam.fetch((&chrom, gene.begin, gene.end)) {
                for read in bam.records() {
                    let asg: bam::Record = match read {
                        Ok(asg) => asg,
                        Err(e) => panic!("Cannot parse bam records {}", e.to_string()),
                    };
                    if asg.is_secondary()
                        || asg.is_supplementary()
                        || asg.is_quality_check_failed()
                        || asg.mapq() < min_q
                        || (asg.is_duplicate() && ignore_duplicate)
                        || gene_strand != asg.is_reverse()
                    {
                        continue;
                    }
                    *gene_type_read_dist
                        .entry(gene.gene_type.to_owned())
                        .or_insert_with(|| HashMap::new())
                        .entry(asg.seq_len() as u16)
                        .or_insert(0) += 1;
                }
            }
        }
    }
    gene_type_read_dist
}

pub fn wrapper_compute_gene_type_read_dist(
    bam_file: &str,
    min_q: u8,
    ignore_duplicate: bool,
    gff3_file: &str,
    features: Vec<String>,
    gene_id: String,
    gene_name: String,
    gene_type: String,
) -> HashMap<String, HashMap<u16, u32>> {
    let gene_map = gff3_parser(gff3_file, features, gene_id, gene_name, gene_type);
    compute_gene_type_read_dist(bam_file, gene_map, min_q, ignore_duplicate)
}
