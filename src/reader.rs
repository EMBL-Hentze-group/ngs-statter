use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{BufReader, Read};
/// Reads a file and returns a buffered reader for iterating over its lines.
///
/// # Arguments
///
/// * `file_path` - A string slice that represents the path to the file.
///
/// # Panics
///
/// This function will panic if the file cannot be found.
///
/// # Returns
///
/// A `BufReader` wrapped in a `Box<dyn Read>`, which can be used to iterate over the lines of the file.
/// If the file is a gzip file (based on its extension), the reader will automatically decode the gzip data.
/// Otherwise, the reader will read the file as is.
///
/// # Examples
///
/// ```
/// use std::io::BufRead;
///
/// let file_path = "path/to/file.txt";
/// let reader = read_file_lines(file_path);
///
/// for line in reader.lines() {
///     println!("{}", line.unwrap());
/// }
/// ```
///
pub fn file_reader(file_path: &str) -> BufReader<Box<dyn Read>> {
    let file: File = match File::open(file_path) {
        Ok(file) => file,
        Err(_) => panic!("Unable to find file {:?}", file_path),
    };

    // Check if the file is a gzip file based on its extension
    let is_gzip: bool = file_path.ends_with(".gz");

    let reader: Box<dyn Read> = if is_gzip {
        Box::new(GzDecoder::new(file))
    } else {
        Box::new(file)
    };
    BufReader::new(reader)
}
