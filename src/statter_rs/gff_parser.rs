use crate::reader;
use regex::Regex;
use std::{
    collections::{HashMap, HashSet},
    io::{BufRead, BufReader, Read},
};

/// Represents a gene extracted from a genomic annotation.
///
/// This struct holds information about a gene, such as its location on a chromosome,
/// start and end positions, strand, unique identifier, name, and type.
#[derive(Debug)]
pub struct Gene {
    /// The chromosome on which the gene is located.
    pub chromosome: String,

    /// The starting position of the gene on the chromosome (0-based).
    pub begin: u32,

    /// The ending position of the gene on the chromosome (exclusive).
    pub end: u32,

    /// The strand on which the gene is located ("+" for forward, "-" for reverse).
    pub strand: String,

    /// The unique identifier of the gene.
    pub gene_id: String,

    /// The name of the gene.
    pub gene_name: String,

    /// The type or category of the gene.
    pub gene_type: String,
}

/// Parses a GFF3 file and extracts information for specified features into a HashMap.
///
/// # Arguments
///
/// * `gff3` - The GFF3 file content as a string.
/// * `features` - A vector of feature types to extract information for.
/// * `gene_id` - The attribute name for gene ID.
/// * `gene_name` - The attribute name for gene name.
/// * `gene_type` - The attribute name for gene type.
///
/// # Returns
///
/// A HashMap where keys are chromosome names, and values are vectors of `Gene` structs.
///
/// # Panics
///
/// Panics if it encounters an error parsing the GFF3 file.
///
/// # Examples
///
/// ```
/// let gff3_content = "your_gff3_content_here";
/// let features = vec!["exon".to_string(), "CDS".to_string()];
/// let gene_id = "ID".to_string();
/// let gene_name = "Name".to_string();
/// let gene_type = "Type".to_string();
///
/// let result = gff3_parser(gff3_content, features, gene_id, gene_name, gene_type);
/// ```
pub fn gff3_parser(
    gff3: &str,
    features: Vec<String>,
    gene_id: String,
    gene_name: String,
    gene_type: String,
) -> HashMap<String, Vec<Gene>> {
    let pattern = Regex::new(r#"(?<k>[^;]+)\=(?<v>[^;]+)"#).unwrap(); // regex
    let mut chrom_map: HashMap<String, Vec<Gene>> = HashMap::new();
    let reader: BufReader<Box<dyn Read>> = reader::file_reader(gff3);
    let feature_set: HashSet<String> = HashSet::from_iter(features);
    let mut feat_counter: u32 = 0;
    for line in reader.lines() {
        let ann = match line {
            Ok(a) => a,
            Err(_) => panic!("Cannot parse the file {:?}", gff3),
        };
        if ann.starts_with("#") {
            continue;
        }
        let annotations: Vec<String> = ann.trim().split("\t").map(|s| s.to_string()).collect();
        if annotations.len() < 9 {
            continue;
        }
        if !feature_set.contains(&annotations[2]) {
            continue;
        }
        let attribs_vec: Vec<(&str, &str)> = pattern
            .captures_iter(&annotations[8])
            .map(|m| {
                let k = m.name("k").unwrap().as_str();
                let v = m.name("v").unwrap().as_str();
                (k, v)
            })
            .collect();
        let mut attribs: HashMap<String, String> = HashMap::new();
        for av in attribs_vec {
            attribs.insert(av.0.to_string(), av.1.to_string());
        }
        let id_gene = match attribs.get(&gene_id) {
            Some(g) => g,
            None => {
                eprintln!(
                    "Cannot find gene_id atrribute {:?} in attribute column {:?}, skipping",
                    gene_id, annotations[8]
                );
                continue;
            }
        };
        let type_gene = match attribs.get(&gene_type) {
            Some(g) => g,
            None => {
                eprintln!(
                    "Cannot find gene_type atrribute {:?} in attribute column {:?}, skipping",
                    gene_type, annotations[8]
                );
                continue;
            }
        };
        let name_gene = match attribs.get(&gene_name) {
            Some(g) => g,
            None => {
                eprintln!(
                    "Cannot find gene_name atrribute {:?} in attribute column {:?}, using 'Unknown'",
                    gene_name, annotations[8]
                );
                "Unknown"
            }
        };
        let begin: u32 = annotations[3].parse::<u32>().unwrap() - 1;
        let end: u32 = annotations[4].parse::<u32>().unwrap();
        let gene = Gene {
            chromosome: annotations[0].to_string(),
            begin: begin,
            end: end,
            strand: annotations[6].to_string(),
            gene_id: id_gene.to_string(),
            gene_type: type_gene.to_string(),
            gene_name: name_gene.to_string(),
        };
        let chr_gene: &mut Vec<Gene> = chrom_map.entry(annotations[0].to_string()).or_default();
        chr_gene.push(gene);
        feat_counter += 1;
    }
    if chrom_map.len() == 0 {
        panic!(
            "Cannot parse {:?} features from gff3 file {}",
            feature_set, gff3
        );
    }
    eprintln!(
        "Parsed {} {:?} features from {:?}",
        feat_counter, feature_set, gff3
    );
    chrom_map
}
