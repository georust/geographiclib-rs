extern crate criterion;
extern crate geographiclib;
extern crate geographiclib_rs;

use criterion::{criterion_group, criterion_main, Criterion};
use std::time::Duration;

use geographiclib_rs::{DirectGeodesic, InverseGeodesic};
use std::fs::File;
use std::io::{BufRead, BufReader};

const TEST_MODE_FULL: &str = "full";
const TEST_MODE_SHORT: &str = "short";
const TEST_MODE_DEFAULT: &str = "default";

const FULL_TEST_PATH: &str = "test_fixtures/test_data_unzipped/GeodTest.dat";
const SHORT_TEST_PATH: &str = "test_fixtures/test_data_unzipped/GeodTest-short.dat";
const BUILTIN_TEST_PATH: &str = "test_fixtures/GeodTest-100.dat";
fn test_input_path() -> (&'static str, &'static str) {
    if cfg!(feature = "test_full") {
        (TEST_MODE_FULL, FULL_TEST_PATH)
    } else if cfg!(feature = "test_short") {
        (TEST_MODE_SHORT, SHORT_TEST_PATH)
    } else {
        (TEST_MODE_DEFAULT, BUILTIN_TEST_PATH)
    }
}

fn geodesic_direct_benchmark(c: &mut Criterion) {
    let (mode, file_path) = test_input_path();
    let file = File::open(file_path).unwrap();
    let reader = BufReader::new(file);
    let inputs: Vec<(f64, f64, f64, f64)> = reader
        .lines()
        .map(|line| {
            let line = line.unwrap();
            let fields: Vec<f64> = line.split(' ').map(|s| s.parse::<f64>().unwrap()).collect();
            (fields[0], fields[1], fields[2], fields[6])
        })
        .collect();

    {
        let mut group = c.benchmark_group("direct (c wrapper)");
        if mode == TEST_MODE_FULL {
            group.measurement_time(Duration::from_secs(30));
        }
        group.bench_function(mode, |b| {
            let geod = geographiclib::Geodesic::wgs84();
            b.iter(|| {
                for (lat1, lon1, azi1, s12) in inputs.clone() {
                    let (_lat2, _lon2, _azi2) = geod.direct(lat1, lon1, azi1, s12);
                }
            })
        });
        group.finish();
    }

    {
        let mut group = c.benchmark_group("direct (rust impl)");
        if mode == TEST_MODE_FULL {
            group.measurement_time(Duration::from_secs(35));
        }
        group.bench_function(mode, |b| {
            let geod = geographiclib_rs::Geodesic::wgs84();
            b.iter(|| {
                for (lat1, lon1, azi1, s12) in inputs.clone() {
                    // Do work comparable to geographiclib c-wrapper's `geod.direct` method
                    let (_lat2, _lon2, _azi2) = geod.direct(lat1, lon1, azi1, s12);
                }
            })
        });
        group.finish();
    }
}

fn geodesic_inverse_benchmark(c: &mut Criterion) {
    let (mode, file_path) = test_input_path();
    let file = File::open(file_path).unwrap();
    let reader = BufReader::new(file);
    let inputs: Vec<(f64, f64, f64, f64)> = reader
        .lines()
        .map(|line| {
            let line = line.unwrap();
            let fields: Vec<f64> = line.split(' ').map(|s| s.parse::<f64>().unwrap()).collect();
            (fields[0], fields[1], fields[3], fields[4])
        })
        .collect();

    {
        let mut group = c.benchmark_group("inverse (c wrapper)");
        if mode == TEST_MODE_FULL {
            group.measurement_time(Duration::from_secs(50));
        }
        group.bench_function(mode, |b| {
            let geod = geographiclib::Geodesic::wgs84();
            b.iter(|| {
                for (lat1, lon1, lat2, lon2) in inputs.clone() {
                    let (_s12, _azi1, _azi2, _a12) = geod.inverse(lat1, lon1, lat2, lon2);
                }
            })
        });
        group.finish();
    }

    {
        let mut group = c.benchmark_group("inverse (rust impl)");
        if mode == TEST_MODE_FULL {
            group.measurement_time(Duration::from_secs(70));
        }
        group.bench_function(mode, |b| {
            let geod = geographiclib_rs::Geodesic::wgs84();
            b.iter(|| {
                for (lat1, lon1, lat2, lon2) in inputs.clone() {
                    // Do work comparable to geographiclib c-wrapper's `geod.inverse` method
                    let (_s12, _azi1, _azi2, _a12) = geod.inverse(lat1, lon1, lat2, lon2);
                }
            })
        });
        group.finish();
    }
}

criterion_group!(
    benches,
    geodesic_direct_benchmark,
    geodesic_inverse_benchmark
);
criterion_main!(benches);
