use rust_htslib::{bam, bam::Read};
use std::collections::HashMap;

/// This code snippet represents a function called `align_stats` that takes two parameters: `bam` (a string) and `min_q` (an unsigned 8-bit integer).
/// It returns a `HashMap` containing statistics about alignment.
/// The function initializes an empty `HashMap` called `map_stat` to store the statistics. It then attempts to create a `bam::Reader` from the given `bam` path.
/// If successful, it iterates over the records in the BAM file.
/// For each record, it checks if it is unmapped using the `is_unmapped` method.
/// If true, it increments the count for the "Unmapped" key in `map_stat`.
/// If mapped, it checks for other conditions such as being secondary, supplementary, quality check failed, or having a mapping quality below `min_q`.
/// If any of these conditions are met, the loop continues to the next record. Otherwise, it increments the count for the "Mapped" key in `map_stat`.
/// Finally, it increments the count for the "Input reads" key in `map_stat` for every record.
/// The function returns the `map_stat` `HashMap` containing the alignment statistics.
/// Example usage:
/// ```
/// let bam_path = "path/to/bam";
/// let min_quality = 20;
/// let stats = align_stats(bam_path, min_quality);
/// println!("{:?}", stats);
/// ```
/// This will print the alignment statistics for the given BAM file.
///
pub fn alignment_stats(bam: &str, min_q: u8) -> HashMap<String, u32> {
    let mut map_stat: HashMap<String, u32> = HashMap::new();
    let mut aln_bam: bam::Reader = match bam::Reader::from_path(bam) {
        Ok(sbam) => sbam,
        Err(e) => panic!("Cannot read {}: {}", bam, e.to_string()),
    };
    for aln in aln_bam.records() {
        let asg: bam::Record = match aln {
            Ok(asg) => asg,
            Err(e) => panic!("Cannot parse bam records {}", e.to_string()),
        };
        if asg.is_secondary()
            || asg.is_supplementary()
            || asg.is_quality_check_failed()
            || asg.mapq() < min_q
        {
            continue;
        } else if asg.is_unmapped() {
            *map_stat.entry(String::from("unmapped")).or_insert(0) += 1;
        } else {
            *map_stat.entry(String::from("mapped")).or_insert(0) += 1;
        }
        *map_stat.entry(String::from("input_reads")).or_insert(0) += 1;
    }
    map_stat
}
