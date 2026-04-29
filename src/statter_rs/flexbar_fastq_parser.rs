use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

/// Fix UMI positions in FASTQ headers produced by Flexbar.
/// Moves the UMI from the comment section to the end of the read identifier.
pub fn fix_umi_header<P: AsRef<Path>>(in_fq: P, out_fq: P, separator: &str) -> io::Result<()> {
    let in_path = in_fq.as_ref();
    let out_path = out_fq.as_ref();
    let reader: Box<dyn BufRead> = if in_path.extension().map_or(false, |e| e == "gz" || e == "gzip") {
        Box::new(BufReader::new(GzDecoder::new(File::open(in_path)?)))
    } else {
        Box::new(BufReader::new(File::open(in_path)?))
    };
    let writer: Box<dyn Write> = if out_path.extension().map_or(false, |e| e == "gz" || e == "gzip") {
        Box::new(BufWriter::new(GzEncoder::new(File::create(out_path)?, Compression::default())))
    } else {
        Box::new(BufWriter::new(File::create(out_path)?))
    };
    let mut writer = BufWriter::new(writer);
    let mut lines = reader.lines();
    while let (Some(Ok(header)), Some(Ok(seq)), Some(Ok(plus)), Some(Ok(qual))) =
        (lines.next(), lines.next(), lines.next(), lines.next())
    {
        if !header.starts_with('@') {
            writeln!(writer, "{}", header)?;
            writeln!(writer, "{}", seq)?;
            writeln!(writer, "{}", plus)?;
            writeln!(writer, "{}", qual)?;
            continue;
        }
        let mut parts = header[1..].splitn(2, ' ');
        let name = parts.next().unwrap_or("");
        let comment = parts.next();
        if let Some(comment) = comment {
            let dat: Vec<&str> = comment.split(separator).collect();
            if dat.len() < 2 {
                writeln!(writer, "{}", header)?;
            } else {
                let umi = dat.last().unwrap();
                let new_header = format!(@"{}{}{} {}", name, separator, umi, dat[..dat.len()-1].join(" "));
                writeln!(writer, "@{}", new_header)?;
            }
        } else {
            writeln!(writer, "@{}", name)?;
        }
        writeln!(writer, "{}", seq)?;
        writeln!(writer, "{}", plus)?;
        writeln!(writer, "{}", qual)?;
    }
    Ok(())
}
