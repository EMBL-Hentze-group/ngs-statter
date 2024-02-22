use rust_htslib::{bam, bam::record::Aux, bam::Read};
use std::collections::HashMap;
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

// return STAR aligner alignment stats as a hashmap
pub fn bam_stats(bam: &str, min_q: u8) -> HashMap<String, u32>{
    let mut map_stat: HashMap<String, u32> = HashMap::new();
    let mut star_bam: bam::Reader = match bam::Reader::from_path(bam) {
        Ok(sbam) => sbam,
        Err(e) => panic!("Cannot read {}: {}", bam, e.to_string()),
    };
    for aln in star_bam.records(){
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
    }
    map_stat
}

