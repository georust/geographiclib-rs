extern crate geographiclib_rs;

use geographiclib_rs::Geodesic;
use std::collections::{HashMap};
use std::error::Error;
use std::fs::{self, File};
use std::io::{self, BufRead, prelude::*};
use std::path;
use std::time;
use std::u64;
use zip;

const ZIP_PATH_RELATIVE: &str = "test_fixtures/geographiclib-instrumented/geographiclib-instrumented-crude-out";
pub const DAT_PATH_RELATIVE: &str = "test_fixtures/test-data-unzipped";

// Expected paths for zipped and unzipped versions of a file.
pub struct DataPathPair {
    pub path_zip: path::PathBuf,
    pub path_dat: path::PathBuf,
}

// Given an operation name, return the path to the zipped and unzipped dat file for that operation.
// Note that the unzipped path may point to a location that doesn't exist.
pub fn get_data_paths(op_name: &str) -> std::io::Result<DataPathPair> {
    let mut filename_zip = op_name.to_owned();
    let mut filename_dat = op_name.to_owned();
    filename_zip.push_str(".dat.zip");
    filename_dat.push_str(".dat");

    let dir_base = std::env::current_dir()?;
    let path_base = dir_base.as_path();
    let path_zip = path::Path::new(path_base).join(ZIP_PATH_RELATIVE);
    let path_dat = path::Path::new(path_base).join(DAT_PATH_RELATIVE);
    let result = DataPathPair {
        path_zip: path_zip.join(filename_zip),
        path_dat: path_dat.join(filename_dat),
    };
    Ok(result)
}

// Given expected paths to the zipped and unzipped copies of a data file,
// check whether the unzipped copy is missing, or whether the zipped copy
// is newer. In either of those cases, unzip a fresh copy of the zipped file.
pub fn refresh_unzipped(paths: &DataPathPair) -> std::io::Result<()> {
    let path_zip = paths.path_zip.as_path();
    let modified_in = path_zip.metadata().expect("Failed to read zip file metadata")
        .modified().expect("Failed to read zip file modified date");
    let path_dat = paths.path_dat.as_path();
    let modified_out = match path_dat.metadata() {
        Ok(meta) => meta.modified().expect("Failed to read dat file metadata"),
        Err(_error) => {
            fs::create_dir_all(path_dat.parent().expect("Failed to create dat folder"))?;
            time::SystemTime::UNIX_EPOCH
        }
    };
    if modified_out < modified_in {
        // The unzipped version either doesn't exist, or is older than the zipped version.
        let file_in = fs::File::open(path_zip).expect("Failed to open zip file");
        let mut zip_archive = zip::ZipArchive::new(file_in).expect("Failed to read zip archive");
        if zip_archive.len() != 1 {
            panic!("Expected zip archive to contain exactly one file");
        }
        let mut zip_first = zip_archive.by_index(0).expect("Failed to open first file in zip archive");
        let mut zip_first_contents = String::new();
        zip_first.read_to_string(&mut zip_first_contents).expect("Failed to read zip file content");
        fs::write(path_dat, zip_first_contents).expect("Failed to write zip content to dat file");
    }
    Ok(())
}

// Convert a string to f64, handling various special cases
// found in geographiclib C++ instrumented-crude data files.
pub fn as_f64(s: &str) -> Result<f64, Box<dyn Error + Sync + Send>> {
    match s {
        "nan" => Ok(f64::NAN),
        "-nan" => Ok(-f64::NAN),
        "inf" => Ok(f64::INFINITY),
        "-inf" => Ok(f64::NEG_INFINITY),
        _ => {
            if s.starts_with("0x") {
                // Bit fields are recorded using hex notation, to avoid endian difference complications.
                let without_prefix = s.trim_start_matches("0x");
                match u64::from_str_radix(without_prefix, 16) {
                    Ok(num) => Ok(num as f64),
                    Err(err) => Err(Box::<dyn Error + Sync + Send>::from(err.to_string())),
                }
            } else {
                // All other value types are recorded using values that can be parsed as f64.
                // In some cases, they'll need to be coerced into a different type later.
                match s.parse::<f64>() {
                    Ok(num) => Ok(num),
                    Err(err) => Err(Box::<dyn Error + Sync + Send>::from(err.to_string())),
                }
            }
        },
    }
}

// Return the absolute difference between two values.
// If both values are nan, consider the difference to be 0.
pub fn calc_delta_abs(x: f64, y: f64) -> f64 {
    if x.is_nan() && y.is_nan() {
        return 0f64;
    } else if x.is_infinite() && y.is_infinite() {
        return if x.is_sign_negative() == y.is_sign_negative() { 0f64 } else { f64::INFINITY };
    } else {
        let diff_abs = (x - y).abs();
        return diff_abs;
    }
}

// Return the absolute and relative difference between two values, in that order.
// If both values are nan, consider the difference to be 0.
pub fn calc_delta(x: f64, y: f64) -> (f64, f64) {
    let diff_abs = calc_delta_abs(x, y);
    let diff_rel = if diff_abs == 0.0 || diff_abs.is_nan() || diff_abs.is_infinite() {
        diff_abs
    } else {
        2.0 * diff_abs / (x.abs() + y.abs())
    };
    (diff_abs, diff_rel)
}

// Return true if delta a is "worse" than delta b.
// NAN is worse than INFINITY is worse than anything finite.
// All deltas are assumed positive.
pub fn is_delta_worse(a: f64, b: f64) -> bool {
    (a.is_nan() && !b.is_nan()) || a > b
}

// Extract a the first value from each line, and return as a vector.
pub fn as_vec_basic(op_name: &str) -> Vec<f64> {
    let paths = get_data_paths(op_name).expect("Failed to determine zip and dat paths");
    refresh_unzipped(&paths).expect("Failed to refresh unzipped data file");
    let file = File::open(paths.path_dat.as_path()).unwrap();
    let reader = io::BufReader::new(file);
    let result: Vec<f64> = reader.lines().enumerate()
        // Skip the header line
        .filter(|(i, _line)| *i != 0)
        .map(|(_i, line)| {
            let line_safe = line.expect("Failed to read line");
            let item = line_safe.splitn(2, ' ').nth(0).unwrap();
            as_f64(item).expect("Failed to parse item")
        })
        .collect();
        result
}

// Construct a lookup of geodesics by construction parameters.
// We can't use f64 in the key directly because it implements PartialEq,
// but not Eq or Hash, because of the complication that NaN != NaN.
// See confirm_geodesic and get_geodesic for other helper functionality
// related to this same type of lookup.
pub fn get_geodesic_lookup(data: &Vec<Vec<f64>>) -> HashMap<(u64, u64), Geodesic> {
    let mut map: HashMap<(u64, u64), Geodesic> = HashMap::new();
    for items in data {
        confirm_geodesic(items[0], items[1], &mut map);
    }
    map
}
// Check whether a geodesic with the given parameters exists in the lookup.
// If none is found, add one.
fn confirm_geodesic(a: f64, f: f64, map: &mut HashMap<(u64, u64), Geodesic>) {
    assert!(!a.is_nan() && !f.is_nan(), "This operation is not suitable for use with NaNs");
    let key: (u64, u64) = (a.to_bits(), f.to_bits());
    if !map.contains_key(&key) {
        map.insert(key, Geodesic::new(a, f));
    }
    // map.get(&key).expect("Failed to retrieve Geodesic")
}
// For use with get_geodesic_lookup
pub fn get_geodesic(a: f64, f: f64, map: &HashMap<(u64, u64), Geodesic>) -> &Geodesic {
    assert!(!a.is_nan() && !f.is_nan(), "This operation is not suitable for use with NaNs");
    let key: (u64, u64) = (a.to_bits(), f.to_bits());
    map.get(&key).expect("Failed to find Geodesic")
}

// Extract a specified number of values from each line,
// and return as a vector of vectors.
pub fn as_vecs_basic(op_name: &str, arg_count: isize) -> Vec<Vec<f64>> {
    let paths = get_data_paths(op_name).expect("Failed to determine zip and dat paths");
    refresh_unzipped(&paths).expect("Failed to refresh unzipped data file");
    let path_dat = paths.path_dat.as_path();
    let file = match File::open(path_dat) {
        Ok(val) => val,
        Err(_error) => {
            let path_str = path_dat.to_str().expect("Failed to data file path to string during error reporting");
            panic!("Failed to open data file {}", path_str)
        }
    };
    let reader = io::BufReader::new(file);
    let result: Vec<Vec<f64>> = reader.lines().enumerate()
        // Skip the header line
        .filter(|(i, _line)| *i != 0)
        .map(|(_i, line)| {
                line.expect("Failed to read line").split(' ')
                .enumerate()
                .filter(|(i, _item)| arg_count < 0 || *i < arg_count as usize)
                .map(|(_i, item)| as_f64(item).expect("Failed to parse item"))
                .collect()
        })
        .collect();
    result
}

// Extract numeric values at specified indices on each line,
// and return as a vector of vectors.
pub fn as_vecs_num_sparse(file: &File, header_count: usize, indices: &[usize]) -> Vec<Vec<f64>> {
    let reader = io::BufReader::new(file);
    let result: Vec<Vec<f64>> = reader.lines().enumerate()
        // Skip header lines
        .filter(|(i, _line)| *i >= header_count)
        .map(|(_i, line)| {
                line.expect("Failed to read line").split(' ')
                .enumerate()
                .filter(|(i, _item)| indices.contains(i))
                .map(|(_i, item)| as_f64(item).expect("Failed to parse item"))
                .collect()
        })
        .collect();
    result
}

// Read an unzipped copy of Karney's full GeodTest.dat file from a standard relative path.
pub fn read_geodtest() -> File {
    let dir_base = std::env::current_dir().expect("Failed to determine current directory");
    let path_base = dir_base.as_path();
    let pathbuf = path::Path::new(path_base)
        .join(DAT_PATH_RELATIVE)
        .join("GeodTest.dat");
    let path = pathbuf.as_path();
    let file = match File::open(path) {
        Ok(val) => val,
        Err(_error) => {
            let path_str = path.to_str().expect("Failed to convert GeodTest path to string during error reporting");
            panic!("Failed to open GeodTest.dat file. It may need to be downloaded and unzipped to: {}", path_str)
        }
    };
    file
}

// Centralized logic for reading a geographiclib instrumented-crude *_consts data file.
// Each such data file should contain a single header line, and a single data line.
// op_name: The base name of the data file to look for.
// arg_count: The number of values expected on the data line, or negative to skip this check.
// Returns a vector of values reflecting the items on the data line.
pub fn read_consts_basic(op_name: &str, arg_count: isize) -> Vec<f64> {
    let paths = get_data_paths(op_name).expect("Failed to determine zip and dat paths");
    refresh_unzipped(&paths).expect("Failed to refresh unzipped data file");
    let file = File::open(paths.path_dat.as_path()).unwrap();
    let reader = io::BufReader::new(file);
    let mut result: Vec<f64> = Vec::new();
    reader.lines().enumerate()
        // Skip the header line
        .filter(|(i, _line)| *i != 0)
        .for_each(|(i, line)| {
            assert!(i == 1, "Expected exactly one data line in consts file, but found multiple.");
            let line_safe = line.expect("Failed to read line");
            result = line_safe.split(' ').enumerate()
                .map(|(j, item)| {
                    match as_f64(item) {
                        Ok(parsed) => parsed,
                        Err(_error) => panic!("Error parsing item {} on line {}: {}", j+1, i+1, item),
                    }                    
                })
                .collect();
            assert!(arg_count < 0 || result.len() == arg_count as usize, "Expected {} items per line. Line {} had {}: {}", arg_count, i, result.len(), line_safe);
        });
        result
}

// Centralized logic for reading through a generic geographiclib instrumented-crude data file.
// Each file has a single header line that indicates operation type, version numbers,
// and data line value meanings. Each additional line is a space-separated list of values.
// op_name: The base name of the data file to look for.
// arg_count: The number of values expected on each data line, or negative to skip this check.
// Calls function f for each data line in the file.
pub fn test_basic<T>(op_name: &str, arg_count: isize, f: T)
 where T: Fn(usize, &Vec<f64>)
{
    let paths = get_data_paths(op_name).expect("Failed to determine zip and dat paths");
    refresh_unzipped(&paths).expect("Failed to refresh unzipped data file");
    let file = File::open(paths.path_dat.as_path()).unwrap();
    test_numeric(&file, 1, arg_count, f);
}

pub fn test_numeric<T>(file: &File, skip_count: usize, arg_count: isize, f: T)
 where T: Fn(usize, &Vec<f64>)
{
    let reader = io::BufReader::new(file);
    reader.lines().enumerate()
        // Skip header lines
        .filter(|(i, _line)| *i >= skip_count)
        .for_each(|(i, line)| {
            let line_safe = line.expect("Failed to read line");
            let items: Vec<f64> = line_safe.split(' ').enumerate()
            .map(|(j, item)| {
                match as_f64(item) {
                    Ok(parsed) => parsed,
                    Err(_error) => panic!("Error parsing item {} on line {}: {}", j+1, i+1, item),
                }                    
            })
            .collect();
        assert!(arg_count < 0 || items.len() == arg_count as usize, "Expected {} items per line. Line {} had {}: {}", arg_count, i, items.len(), line_safe);
            // Report 1-based line number, rather than 0-based
            f(i+1, &items);
        });
}

// When displaying f64, stable Rust declines to display the sign for -0 or -nan,
// as of 2021/01/10. It looks like a fix for this is on the way:
//     https://github.com/rust-lang/rust/issues/20596
// Here's a work-around in the meantime. It's mechanics aren't great, but it's temporary.
pub fn help_sign(x: f64) -> String {
    if (x == 0.0 || x.is_nan()) && x.is_sign_negative() {
        "-".to_string()
    } else {
        "".to_string()
    }
}
