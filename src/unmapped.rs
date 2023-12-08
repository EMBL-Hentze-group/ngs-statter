use flate2::{write, Compression};
use rust_htslib::{bam, bam::record::Aux, bam::Read};
use std::collections::HashMap;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

// fastq writer either plain file or gz compressed
// Use the presence of .gz suffix to decide
// source: https://users.rust-lang.org/t/write-to-normal-or-gzip-file-transparently/35561/2
fn fq_writer(fq: &str) -> Box<dyn Write> {
    let fq_path = Path::new(fq);
    let file = match File::create(&fq_path) {
        Ok(f) => f,
        Err(e) => panic!("Cannot create {}: {}", fq_path.display(), e.to_string()),
    };
    if fq_path.extension() == Some(OsStr::new("gz")) {
        Box::new(BufWriter::with_capacity(
            128 * 1024,
            write::GzEncoder::new(file, Compression::default()),
        ))
    } else {
        Box::new(BufWriter::with_capacity(128 * 1024, file))
    }
}

// ** STAR Unmapped reads **
// uT : for unmapped reads, reason for not mapping:
//  0 : no acceptable seed/windows, ”Unmapped other” in the Log.final.out
//  1 : best alignment shorter than min allowed mapped length, ”Unmapped: too short” in the Log.final.out
//  2 : best alignment has more mismatches than max allowed number of mismatches, ”Un-mapped: too many mismatches” in the Log.final.out
//  3 : read maps to more loci than the max number of multimappng loci, ”Multimapping: mapped to too many loci” in the Log.final.out
//  4 : unmapped mate of a mapped paired-end read
//
// write "3" into one file and write all 0-2 into another file

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

// parse bam file from single end data
// split unmapped reads into two fastq files:  multimappers and all others
pub fn unmapped_writer_single_end(bam: &str, fq_multi: &str, fq_other: &str) -> HashMap<String, u32> {
    let mut bam = match bam::Reader::from_path(bam) {
        Ok(bam) => bam,
        Err(e) => panic!("Cannot read {}:{}", bam, e.to_string()),
    };
    let mut multi_writer = fq_writer(fq_multi);
    let mut other_writer = fq_writer(fq_other);
    let mut unmapped_stats: HashMap<String, u32> = HashMap::new();
    for aln in bam.records() {
        let rc = match aln {
            Ok(rc) => rc,
            Err(e) => panic!("Problem parsing bam records {}", e.to_string()),
        };
        if !rc.is_unmapped() {
            continue;
        }
        let ut = match rc.aux(b"uT") {
            // see documentation above for this flag
            Ok(v) => match v {
                Aux::Char(v) => v as char,
                _ => panic!("Expected Aux::Char"),
            },
            Err(_) => continue,
        };
        let mut fq: Vec<u8> = Vec::new();
        fq.push(64); // @
        fq.extend(rc.qname()); // read id
        fq.push(10); // \n
        fq.extend(rc.seq().as_bytes()); // read sequence
        fq.extend(vec![10, 43, 10]); // \n+\n
        let qual: Vec<u8> = rc.qual().into_iter().map(|q| q + 33).collect();
        fq.extend(qual); // PHRED score
        fq.push(10); // new line
        if ut == '3' {
            multi_writer.write_all(&fq).unwrap();
        } else {
            other_writer.write_all(&fq).unwrap();
        }
        unmapped_stats
            .entry(unmapped_kind(&ut))
            .and_modify(|e| *e += 1)
            .or_insert(1);
    }
    unmapped_stats
}