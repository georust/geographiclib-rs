extern crate geographiclib_rs;

use geographiclib_rs::Geodesic;
use std::collections::{BTreeMap, HashMap};
use std::error::Error;
use std::fmt::Display;
use std::fs::{self, File};
use std::io::{self, BufRead, prelude::*};
use std::path;
use std::time;
use std::u64;
use zip;

// testing reminder: Use "cargo test -- --nocapture" to let tests print to console.
// testing reminder: If a new "*_vs_cpp" test fails especially badly, review the C++ "instrumented-crude" logging for possible bugs there
// testing tip: While assert variations are valuable for catching regressions, they're
//              often of more limited value when trying to _improve_ result precision.
//              In such cases, use of something like DeltaEntry is often helpful,
//              but remember to use nocapture, and often directing output to a file, for
//              syntax like "cargo test -- --nocapture > ../test-out-a.txt".
// benchmarking reminder: Remember to shut down unneeded processes before benching.
// benchmarking reminder: In some cases, it may be desirable to capture multiple before and after
//                        bench runs, to get a clearer sense of natural variability.

const ZIP_PATH_RELATIVE: &str = "test_fixtures/geographiclib-instrumented/geographiclib-instrumented-crude-out";
pub const DAT_PATH_RELATIVE: &str = "test_fixtures/test-data-unzipped";

#[allow(non_upper_case_globals)]
pub const nC_: usize = 7; // todo: define in geodesic.rs instead

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

// Return the absolute and relative difference between two values, in that order.
// If both values are nan, consider the difference to be 0.
pub fn calc_delta(x: f64, y: f64) -> (f64, f64) {
    if x.is_nan() && y.is_nan() {
        return (0f64, 0f64);
    } else if x.is_infinite() && y.is_infinite() {
        return if x.is_sign_negative() == y.is_sign_negative() { (0f64, 0f64) } else { (f64::INFINITY, f64::INFINITY) };
    } else {
        let diff_abs = (x - y).abs();
        let diff_rel = if diff_abs == 0.0 {
            0.0
        } else {
            2.0 * diff_abs / (x.abs() + y.abs())
        };
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

// A struct for taking difference values, splitting into buckets, and reporting back.
pub struct DeltaHistogram {
    pub num_nan: usize,
    pub num_inf: usize,
    pub num_zero: usize,
    pub exp_buckets: HashMap<isize, usize>,
}

impl DeltaHistogram {
    pub fn new() -> Self {
        DeltaHistogram {
            num_nan: 0,
            num_inf: 0,
            num_zero: 0,
            exp_buckets: HashMap::new(),
        }
    }

    // Add a new difference item to the dataset being tracked.
    pub fn add(&mut self, delta: f64) {
        if delta.is_nan() {
            self.num_nan += 1;
        } else if delta.is_infinite() {
            self.num_inf += 1;
        } else if delta == 0.0 {
            self.num_zero += 1;
        } else {
            let exp = delta.log10() as isize;
            let current: usize = match self.exp_buckets.get(&exp) {
                Some(val) => *val,
                _ => 0,
            };
            self.exp_buckets.insert(exp, current + 1);
        }
    }
}

impl Clone for DeltaHistogram {
    fn clone(&self) -> Self {
        DeltaHistogram {
            num_nan: self.num_nan,
            num_inf: self.num_inf,
            num_zero: self.num_zero,
            exp_buckets: self.exp_buckets.clone(),
        }
    }
}

// Round a value for use in DeltaHistogram display.
// Never round to 0 or 100. Only accept those values naturally.
fn to_percent(num_part: usize, num_all: usize) -> usize {
    let percent = 100f64 * num_part as f64 / num_all as f64;
    let rounded = if percent < 1.0 && num_part != 0 {
        1
    } else if percent > 99.0 && num_part != num_all {
        99
    } else {
        percent.round() as usize
    };
    rounded
}

impl Display for DeltaHistogram {
    // Display a summary, reduced down to a manageable number of buckets.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        const MAX_BUCKETS: usize = 7;
        let mut keys_asc: Vec<isize> = Vec::new();
        let mut collapsed_from: Vec<isize> = Vec::new();
        let mut collapsed_to: Vec<isize> = Vec::new();
        let mut histo_reduced: BTreeMap<isize, usize> = BTreeMap::new();
        let mut num_total = 0;
        self.exp_buckets.iter().for_each(|(&key, &val)| {
            keys_asc.push(key);
            histo_reduced.insert(key, val);
            num_total += val;
        });
        keys_asc.sort();
        while histo_reduced.len() > MAX_BUCKETS {
            // Collapse the smallest bucket into its least-populated neighbor.
            let mut smallest = usize::MAX;
            let mut collapse_from = isize::MIN;
            self.exp_buckets.iter().for_each(|(&key, &val)| {
                if val < smallest {
                    collapse_from = key;
                    smallest = val;
                }
            });
            let index_smallest = keys_asc.iter().position(|&val| val == collapse_from).unwrap();
            let key_prev = if index_smallest == 0 { isize::MIN } else { keys_asc[index_smallest - 1] };
            let key_next = if index_smallest >= histo_reduced.len() - 1 { isize::MIN } else { keys_asc[index_smallest + 1] };
            let size_prev: usize = if key_prev == isize::MIN {
                usize::MAX
            } else {
                *histo_reduced.get(&key_prev).unwrap()
            };
            let size_next = if key_prev == isize::MIN {
                usize::MAX
            } else {
                *histo_reduced.get(&key_next).unwrap()
            };
            let (collapse_to, size_to) = if size_next < size_prev {
                (key_next, size_next)
            } else {
                (key_prev, size_prev)
            };
            histo_reduced.insert(collapse_to, size_to + smallest);
            histo_reduced.remove(&collapse_from);
            keys_asc.remove(index_smallest);
            collapsed_from.push(collapse_from);
            collapsed_to.push(collapse_to);
        }

        write!(f, "count {}", num_total)?;

        if self.num_zero > 0 {
            let percent_zero = to_percent(self.num_zero, num_total); 
            write!(f, ", zero {}%", percent_zero)?;
        }

        // Convert counts to percentages
        histo_reduced.iter_mut().for_each(|(_key, val)| {
            assert!(*val != 0, "Internal error: Bucket contains no items");
            *val = to_percent(*val, num_total);
        });
        for (key, val) in &histo_reduced {
            let gained = collapsed_to.contains(&key);
            write!(f, ", {}e{} {}%", if gained { "~" } else { "" }, key, val)?;
        }
        if self.num_inf > 0 {
            let percent_inf = to_percent(self.num_inf, num_total);
            write!(f, ", inf {}%", percent_inf)?;
        }
        if self.num_nan > 0 {
            let percent_nan = to_percent(self.num_nan, num_total);
            write!(f, ", nan {}%", percent_nan)?;
        }
        Ok(())
    }
}


// Information about value differences.
// This is typically used to record the worst-case differences found
// among sets of comparable expected vs found values.
pub struct DeltaEntry {
    pub diff_abs: f64,
    pub diff_rel: f64,
    pub num_sign_change_only: usize,
    pub num_sign_change_plus_diff: usize,
    pub line_abs: usize,
    pub line_rel: usize,
    pub histo: DeltaHistogram,
}

impl DeltaEntry {
    pub fn new() -> Self {
        DeltaEntry {
            diff_abs: 0.0,
            diff_rel: 0.0,
            num_sign_change_only: 0,
            num_sign_change_plus_diff: 0,
            line_abs: 0,
            line_rel: 0,
            histo: DeltaHistogram::new(),
        }
    }

    // Given x and y, calculate their absolute difference, relative difference,
    // and sign change status, then check whether any of those values is the
    // worst seen so far for comparable operations. If it is, record the line
    // number and the new worst difference.
    // Note that use of this function typically suggests that the values in question
    // are suspected to not be an ideal match for the expected values.
    pub fn add(&mut self, x: f64, y: f64, line_num: usize) {
        let delta = calc_delta(x, y);
        if (f64::is_nan(delta.0) && !f64::is_nan(self.diff_abs)) || delta.0 > self.diff_abs {
            self.diff_abs = delta.0;
            self.line_abs = line_num;
        }
        if (f64::is_nan(delta.1) && !f64::is_nan(self.diff_rel)) || delta.1 > self.diff_rel {
            self.diff_rel = delta.1;
            self.line_rel = line_num;
        }
        if x.is_sign_negative() == y.is_sign_negative() {
        } else if delta.0 == 0.0 {
            self.num_sign_change_only += 1;
        } else {
            self.num_sign_change_plus_diff += 1;
        }
        // On the next line, negating the inequality is intentional, to get desired nan behavior.
        let diff_better = if delta.0.is_nan() || !(delta.0 <= delta.1) {
            delta.1
        } else {
            delta.0
        };
        self.histo.add(diff_better);
    }
}

impl Clone for DeltaEntry {
        fn clone(&self) -> Self {
            DeltaEntry {
                diff_abs: self.diff_abs,
                diff_rel: self.diff_rel,
                num_sign_change_only: self.num_sign_change_only,
                num_sign_change_plus_diff: self.num_sign_change_plus_diff,
                line_abs: self.line_rel,
                line_rel: self.line_rel,
                histo: self.histo.clone(),
            }
        }
}

impl Display for DeltaEntry {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(
            f,
            "max abs {:e} at line {}, max rel {:e} at line {}, sign diffs {}+{}, histo {}",
            self.diff_abs, self.line_abs,
            self.diff_rel, self.line_rel,
            self.num_sign_change_only, self.num_sign_change_plus_diff,
            self.histo,
        )?;
        Ok(())
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
