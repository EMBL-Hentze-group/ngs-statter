use pyo3::panic;
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

/// This function `bam_stats` takes a BAM file path (`bam`) and a minimum quality score (`min_q`) as input parameters.
/// It returns a `HashMap` containing statistics about the alignments in the BAM file.
/// The function reads the BAM file using the `bam::Reader` and iterates over each alignment record.
/// If the alignment is unmapped, it extracts the `uT` tag value and determines the kind of unmapped alignment using the `unmapped_kind` function.
/// It then increments the corresponding entry in the `map_stat` HashMap.
/// If the alignment is mapped, it checks for various conditions such as secondary alignment, supplementary alignment, failed quality check, and minimum mapping quality score.
/// If any of these conditions are met, the function continues to the next alignment record.
/// If the alignment is a PCR duplicate, it increments the corresponding entry in the `map_stat` HashMap.
/// Otherwise, it increments the entry for unique alignments.
/// The function also extracts the `NH` tag value to determine if the alignment is multimapped or unique.
/// Based on this information, it increments the corresponding entry in the `map_stat` HashMap.
/// Finally, the function increments the entries for mapped alignments and input reads.
/// The `map_stat` HashMap is then returned as the result.
///
pub fn bam_stats(bam: &str, min_q: u8) -> HashMap<String, u32> {
    let mut map_stat: HashMap<String, u32> = HashMap::new();
    let mut star_bam: bam::Reader = match bam::Reader::from_path(bam) {
        Ok(sbam) => sbam,
        Err(e) => panic!("Cannot read {}: {}", bam, e.to_string()),
    };
    map_stat.insert(String::from("Mapped: PCR duplicate reads"), 0);
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
            *map_stat.entry(String::from("Unmapped: Total")).or_insert(0) += 1;
            *map_stat.entry(String::from("Input reads")).or_insert(0) += 1;
            continue;
        } else if asg.is_secondary()
            || asg.is_supplementary()
            || asg.is_quality_check_failed()
            || asg.mapq() < min_q
        {
            continue;
        } else if asg.is_duplicate() {
            *map_stat
                .entry(String::from("Mapped: PCR duplicate reads"))
                .or_insert(0) += 1;
        } else {
            *map_stat
                .entry(String::from("Mapped: Unique reads"))
                .or_insert(0) += 1;
        }
        let nh = match asg.aux(b"NH") {
            Ok(v) => match v {
                Aux::U8(v) => v as u32,
                _ => panic!("Expected Aux::U8"),
            },
            Err(_) => panic!("Cannot parse alignment info!"),
        };
        if nh > 1 {
            *map_stat
                .entry(String::from("Mapped: Multimapped reads"))
                .or_insert(0) += 1;
        } else {
            *map_stat
                .entry(String::from("Mapped: Uniquely mapped reads"))
                .or_insert(0) += 1;
        }
        *map_stat.entry(String::from("Mapped: Total")).or_insert(0) += 1;
        *map_stat.entry(String::from("Input reads")).or_insert(0) += 1;
    }
    map_stat
}
