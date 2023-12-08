use pyo3::prelude::*;
use std::collections::HashMap;
mod unmapped;

/// parse single end bam file and
/// split unmapped fastqs into two files: (i) all multimappers (ii) other unmapped reads
/// generate unmapped statistics as json
#[pyfunction]
fn single_end_unmapped_fastq(bam: &str, fq_multi: &str, fq_other: &str) -> HashMap<String, u32>{
    let unmapped_stat = unmapped::unmapped_writer_single_end(bam, fq_multi, fq_other);
    unmapped_stat
}
/// A Python module implemented in Rust.
#[pymodule]
fn statter(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(single_end_unmapped_fastq, m)?)?;
    Ok(())
}
