use rust_htslib::{bam, bam::record::Aux, bam::Read};
use std::collections::HashMap;
// use std::collections::{HashMap, HashSet};
// use std::fs::File;
// use std::io::{BufRead, BufReader};
// read bam file created using STAR and parse alignment stats

// return unmapped kind, source: STAR documentation
// see https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf section Unmapped reads
fn unmapped_kind(ut: &char) -> String {
    match ut {
        '0' => String::from("Unmapped: no seed/windows"),
        '1' => String::from("Unmapped: too short"),
        '2' => String::from("Unmapped: too many mismatches"),
        '3' => String::from("Multimapping: mapped to too many loci"),
        '4' => String::from("Unmapped: paired-end mate"),
        _ => String::from("Unknown"),
    }
}

// read chromosomes

// Given a plain text file with one chromosome per line or
// a fasta index file, parse chromosome names and return
// fn chromosome_names(txt: &str) -> HashSet<String> {
//     let file = match File::open(txt) {
//         Ok(file) => file,
//         Err(_) => panic!("Unable to find {:?}", txt),
//     };
//     let reader: BufReader<File> = BufReader::new(file);
//     let mut chroms: HashSet<String> = HashSet::new();
//     for line in reader.lines() {
//         let dat: Vec<String> = line
//             .unwrap()
//             .trim()
//             .split("\t")
//             .map(|s| s.to_string())
//             .collect();
//         chroms.insert(dat[0].clone());
//     }
//     chroms
// }

// return STAR aligner alignment stats as a hashmap
pub fn bam_stats(bam: &str, min_q: u8) -> HashMap<String, u32> {
    let mut map_stat: HashMap<String, u32> = HashMap::new();
    let mut star_bam: bam::Reader = match bam::Reader::from_path(bam) {
        Ok(sbam) => sbam,
        Err(e) => panic!("Cannot read {}: {}", bam, e.to_string()),
    };
    for aln in star_bam.records() {
        let asg: bam::Record = match aln {
            Ok(asg) => asg,
            Err(e) => panic!("Cannot parse bam records {}", e.to_string()),
        };
        if asg.is_unmapped() {
            let ut = match asg.aux(b"uT") {
                Ok(v) => match v {
                    Aux::Char(v) => v as char,
                    _ => panic!("Expected Aux::Char!"),
                },
                Err(_) => continue,
            };
            *map_stat.entry(unmapped_kind(&ut)).or_insert(0) += 1;
            *map_stat.entry(String::from("Input reads")).or_insert(0) += 1;
            continue;
        } else if asg.is_secondary()
            || asg.is_supplementary()
            || asg.is_quality_check_failed()
            || asg.mapq() < min_q
        {
            continue;
        } else if asg.is_duplicate() {
            *map_stat.entry(String::from("PCR duplicate")).or_insert(0) += 1;
        } else {
            *map_stat.entry(String::from("Unique")).or_insert(0) += 1;
        }
        *map_stat.entry(String::from("Total aligned")).or_insert(0) += 1;
        *map_stat.entry(String::from("Input reads")).or_insert(0) += 1;
    }
    map_stat
}
