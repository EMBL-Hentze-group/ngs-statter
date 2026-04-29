use crate::reader;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read};
/// Calculates the length of each read in a FASTQ file and returns a HashMap
/// with the read length as the key and the count as the value.
///
/// # Arguments
///
/// * `fq_file` - A string slice that holds the path to the FASTQ file.
///
/// # Returns
///
/// A HashMap<u16, u32> where the key is the read length and the value is the count
/// of reads with that length.
///
/// # Panics
///
/// This function will panic if it encounters any errors while parsing the FASTQ file.
///
pub fn fastq_read_length_counter(fastq: &str) -> HashMap<u16, u32> {
    let fq_reader: BufReader<Box<dyn Read>> = reader::file_reader(fastq);
    let mut read_stats: HashMap<u16, u32> = HashMap::new();
    let mut fa_index: usize = 0;
    let mut found_header: bool = false;
    for (index, line) in fq_reader.lines().enumerate() {
        // if (index != fa_index) && (index != fa_index - 1) {
        //     // either the sequence line, or the header line (just before sequence)
        //     continue;
        // }
        let seq = match line {
            Ok(s) => s,
            Err(_) => panic!("Cannot parse fastq lines from file {:?}", fastq),
        };
        if seq.starts_with("@") {
            // header line, confirm it starts with "@"
            found_header = true;
            fa_index = index;
            continue;
        }
        if (found_header) && (index == fa_index+1){
            // sequence line and the line before is confirmed to be header
            let slen: u16 = seq.len() as u16;
            fa_index = index + 4;
            *read_stats.entry(slen).or_insert(0) += 1;
            found_header = false;
        }
    }
    read_stats
}
