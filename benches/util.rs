use std::fs::File;
use std::io::{BufRead, BufReader};

pub const TEST_MODE_FULL: &str = "full";
pub const TEST_MODE_SHORT: &str = "short";
pub const TEST_MODE_DEFAULT: &str = "default";

const FULL_TEST_PATH: &str = "test_fixtures/test_data_unzipped/GeodTest.dat";
const SHORT_TEST_PATH: &str = "test_fixtures/test_data_unzipped/GeodTest-short.dat";
const BUILTIN_TEST_PATH: &str = "test_fixtures/GeodTest-100.dat";

#[derive(Clone)]
pub struct TestInput {
    pub lat1: f64,
    pub lon1: f64,
    pub azi1: f64,
    pub lat2: f64,
    pub lon2: f64,
    pub s12: f64,
}

pub fn load_test_data() -> (&'static str, Vec<TestInput>) {
    let (mode, file_path) = if cfg!(feature = "test_full") {
        (TEST_MODE_FULL, FULL_TEST_PATH)
    } else if cfg!(feature = "test_short") {
        (TEST_MODE_SHORT, SHORT_TEST_PATH)
    } else {
        (TEST_MODE_DEFAULT, BUILTIN_TEST_PATH)
    };

    let file = File::open(file_path).unwrap();
    let reader = BufReader::new(file);
    let records = reader
        .lines()
        .map(|line| {
            let line = line.unwrap();
            let f: Vec<f64> = line.split(' ').map(|s| s.parse::<f64>().unwrap()).collect();
            TestInput {
                lat1: f[0],
                lon1: f[1],
                azi1: f[2],
                lat2: f[3],
                lon2: f[4],
                s12: f[6],
            }
        })
        .collect();

    (mode, records)
}
