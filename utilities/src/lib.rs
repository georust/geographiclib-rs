extern crate geographiclib_rs;

use geographiclib_rs::Geodesic;
use std::collections::HashMap;
use std::error::Error;
use std::fs::{self, File};
use std::io::{self, BufRead, prelude::*};
use std::path;
use std::time;
use std::u64;
use zip;

// todo: add a rust unit test that mirrors C++ MGRS::Check

const ZIP_PATH_RELATIVE: &str = "test_fixtures/geographiclib-instrumented/geographiclib-instrumented-crude-out";
const DAT_PATH_RELATIVE: &str = "test_fixtures/test-data-unzipped";

#[allow(non_upper_case_globals)]
pub const nC_: usize = 7; // todo: define in geodesic-rs instead

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

// Return the absolute and relative difference between two values, in that order.
// If both values are nan, consider the difference to be 0.
pub fn calc_delta(x: f64, y: f64) -> (f64, f64) {
    if x.is_nan() && y.is_nan() {
        return (0f64, 0f64);
    } else if x.is_infinite() && y.is_infinite() {
        return if x.is_sign_negative() == y.is_sign_negative() { (0f64, 0f64) } else { (f64::INFINITY, f64::INFINITY) };
    } else {
        let diff_abs = (x - y).abs();
        let diff_rel = 2.0 * diff_abs / (x.abs() + y.abs());
        return (diff_abs, diff_rel);
    }
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

// We can't use f64 in the key directly because it implements PartialEq,
// but not Eq or Hash, because of the complication that NaN != NaN.
pub fn get_geodesic_lookup(data: &Vec<Vec<f64>>) -> HashMap<(u64, u64), Geodesic> {
    let mut map: HashMap<(u64, u64), Geodesic> = HashMap::new();
    for items in data {
        confirm_geodesic(items[0], items[1], &mut map);
    }
    map
}
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
    let file = File::open(paths.path_dat.as_path()).unwrap();
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
    let reader = io::BufReader::new(file);
    reader.lines().enumerate()
        // Skip the header line
        .filter(|(i, _line)| *i != 0)
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

pub struct RecordDeltaEntry {
    pub diff_abs: f64,
    pub diff_rel: f64,
    pub sign_change: bool,
    pub line_abs: usize,
    pub line_rel: usize,
}

impl RecordDeltaEntry {
    pub fn new() -> Self {
        RecordDeltaEntry {
            diff_abs: 0.0,
            diff_rel: 0.0,
            sign_change: false,
            line_abs: 0,
            line_rel: 0,
        }
    }
}

impl Copy for RecordDeltaEntry {}

impl Clone for RecordDeltaEntry {
        fn clone(&self) -> Self {
            RecordDeltaEntry {
                diff_abs: self.diff_abs,
                diff_rel: self.diff_rel,
                sign_change: self.sign_change,
                line_abs: self.line_rel,
                line_rel: self.line_rel,
            }
        }
}

pub fn record_delta(x: f64, y: f64, line_num: usize, entry: &mut RecordDeltaEntry) {
    let delta2 = calc_delta(x, y);
    if (f64::is_nan(delta2.0) && !f64::is_nan(entry.diff_abs)) || delta2.0 > entry.diff_abs {
        entry.diff_abs = delta2.0;
        entry.line_abs = line_num;
    }
    if (f64::is_nan(delta2.1) && !f64::is_nan(entry.diff_rel)) || delta2.1 > entry.diff_rel {
        entry.diff_rel = delta2.1;
        entry.line_rel = line_num;
    }
    if x.is_sign_negative() != y.is_sign_negative() {
        entry.sign_change = true;
    }
}

// Check whether x and y are within delta of each other, considering both absolute and relative difference.
// Consider them close enough if absolute difference is within delta inclusive,
// or if relative difference is within delta exclusive.
// nan values are considered equal to each other, and outside delta from everything else.
#[macro_export]
macro_rules! assert_delta {
    ($x:expr, $y:expr, $delta:expr, $allow_sign_change:expr, $msg:expr, $line_num:expr) => {
            // todo: use the same formatting behavior as "assert!", instead of line_num hack
            // todo: consider using separate deltas for relative and absolute
            let delta2 = utilities::calc_delta($x, $y);
            assert!(
                    delta2.0 <= $delta || delta2.1 < $delta,
                    "assert_delta failed line {}, {}: {:e} vs {:e} diff {:e} rel {:e} outside inclusive {:e}", $line_num, $msg, $x, $y, delta2.0, delta2.1, $delta
            );
            // For the sign change check, allow (NAN vs NAN), but not (0.0 vs -0.0) or (NAN vs -NAN).
            assert!(
                $allow_sign_change || $x.is_sign_negative() == $y.is_sign_negative(),
                "assert_delta failed line {}, {}: {:e} vs {:e} sign difference disallowed", $line_num, $msg, $x, $y
            );
    }
}



// fn unescape_all(s: &String) -> String {
//     let chars1 = s.as_bytes();
//     let mut chars2 = chars1.to_owned();
//     let len = chars1.len();
//     let mut i_from: usize = 0;
//     let mut i_to: usize = 0;
//     const BACKSLASH: u8 = '\\' as u8;
//     const LOWER_S: u8 = 's' as u8;
//     const SPACE: u8 = ' ' as u8;
//     const ZERO: u8 = '0' as u8;
//     const NOTHING: u8 = '\0' as u8;
//     while i_from < len {
//         let c = chars1[i_from];
//         i_from += 1;
//         if c == BACKSLASH {
//             if i_from >= len {
//                 panic!("Invalid escape sequence. Item ended with single backslash.");
//             }
//             let c_swap = match chars1[i_from] {
//                 BACKSLASH => BACKSLASH,
//                 LOWER_S => SPACE,
//                 ZERO => NOTHING,
//                 _ => panic!("Invalid escape sequence \\{}", &s[i_from..i_from+1])
//             };
//             if c_swap != NOTHING {
//                 chars2[i_to] = c_swap;
//                 i_to += 1;
//             }
//         } else {
//             chars2[i_to] = c;
//             i_to += 1;
//         }
//     }

//     for i in i_to..len {
//         chars2[i] = NOTHING;
//     }
//     let copy = String::from_utf8(chars2).expect("Failed to create new string");
//     return copy;
// }
