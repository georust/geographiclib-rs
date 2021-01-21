extern crate criterion;
extern crate geographiclib;
extern crate geographiclib_rs;

use criterion::{criterion_group, criterion_main, Criterion};

use geographiclib_rs::{DirectGeodesic, InverseGeodesic};
use std::fs::File;
use std::io::{BufRead, BufReader};

const FULL_TEST_PATH: &str = &"test_fixtures/test_data_unzipped/GeodTest.dat";
const SHORT_TEST_PATH: &str = "test_fixtures/test_data_unzipped/GeodTest-short.dat";
const BUILTIN_TEST_PATH: &str = "test_fixtures/GeodTest-100.dat";
fn test_input_path() -> &'static str {
    if cfg!(feature = "test_full") {
        FULL_TEST_PATH
    } else if cfg!(feature = "test_short") {
        SHORT_TEST_PATH
    } else {
        BUILTIN_TEST_PATH
    }
}

fn geodesic_direct_benchmark(c: &mut Criterion) {
    let file = File::open(test_input_path()).unwrap();
    let reader = BufReader::new(file);
    let inputs: Vec<(f64, f64, f64, f64)> = reader
        .lines()
        .map(|line| {
            let line = line.unwrap();
            let fields: Vec<f64> = line.split(" ").map(|s| s.parse::<f64>().unwrap()).collect();
            (fields[0], fields[1], fields[2], fields[6])
        })
        .collect();

    c.bench_function("direct (c wrapper)", |b| {
        let geod = geographiclib::Geodesic::wgs84();
        b.iter(|| {
            for (lat1, lon1, azi1, s12) in inputs.clone() {
                let (_lat2, _lon2, _azi2) = geod.direct(lat1, lon1, azi1, s12);
            }
        })
    });

    c.bench_function("direct (rust impl)", |b| {
        let geod = geographiclib_rs::Geodesic::wgs84();
        b.iter(|| {
            for (lat1, lon1, azi1, s12) in inputs.clone() {
                // Do work comparable to geographiclib c-wrapper's `geod.direct` method
                let (_lat2, _lon2, _azi2) = geod.direct(lat1, lon1, azi1, s12);
            }
        })
    });
}

fn geodesic_inverse_benchmark(c: &mut Criterion) {
    let file = File::open(test_input_path()).unwrap();
    let reader = BufReader::new(file);
    let inputs: Vec<(f64, f64, f64, f64)> = reader
        .lines()
        .map(|line| {
            let line = line.unwrap();
            let fields: Vec<f64> = line.split(" ").map(|s| s.parse::<f64>().unwrap()).collect();
            (fields[0], fields[1], fields[3], fields[4])
        })
        .collect();

    c.bench_function("inverse (c wrapper)", |b| {
        let geod = geographiclib::Geodesic::wgs84();
        b.iter(|| {
            for (lat1, lon1, lat2, lon2) in inputs.clone() {
                let (_s12, _azi1, _azi2, _a12) = geod.inverse(lat1, lon1, lat2, lon2);
            }
        })
    });

    c.bench_function("inverse (rust impl)", |b| {
        let geod = geographiclib_rs::Geodesic::wgs84();
        b.iter(|| {
            for (lat1, lon1, lat2, lon2) in inputs.clone() {
                // Do work comparable to geographiclib c-wrapper's `geod.inverse` method
                let (_s12, _azi1, _azi2, _a12) = geod.inverse(lat1, lon1, lat2, lon2);
            }
        })
    });
}

criterion_group!(
    benches,
    geodesic_direct_benchmark,
    geodesic_inverse_benchmark
);
criterion_main!(benches);
