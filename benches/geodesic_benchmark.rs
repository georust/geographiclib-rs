#![allow(non_snake_case)]
#![allow(non_snake_case)]

extern crate criterion;
extern crate geographiclib;
extern crate geographiclib_rs;

use criterion::{criterion_group, criterion_main, Criterion};
use geographiclib_rs::{DirectGeodesic, InverseGeodesic};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::sync::{Arc, Mutex};
use std::time::Duration;
use utilities::{nC_, util};

fn geodesic_direct_benchmark(c: &mut Criterion) {
    const INPUT_FOR_DIRECT: &'static str = "test_fixtures/GeodTest-100.dat";

    let file = File::open(INPUT_FOR_DIRECT).unwrap();
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
    const INPUT_FOR_INVERSE: &'static str = "test_fixtures/GeodTest-100.dat";

    let file = File::open(INPUT_FOR_INVERSE).unwrap();
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

// Benchmarks using Karney's full GeodTest.dat file: *_geodtest_*

fn benchmark_geodtest_geodesic_direct12(c: &mut Criterion) {
    let mut group = c.benchmark_group("Geodesic.direct GeodTest adjusted");
    group.measurement_time(Duration::from_secs_f64(40.0));
    let data_arc = Arc::new(Mutex::new((
        util::as_vecs_num_sparse(&util::read_geodtest(), 0, &[0, 1, 2, 6]),
        geographiclib_rs::Geodesic::wgs84(),
    )));
    group.bench_function("Geodesic.direct (rust impl, GeodTest.dat inputs)", |b| {
        let pair = data_arc.lock().unwrap();
        let data = &pair.0;
        let g = &pair.1;
        b.iter(|| {
            for line in data {
                let (_lat2, _lon2, _azi2, _m12, _M12, _M21, _S12, _a12) =
                    g.direct(line[0], line[1], line[2], line[3]);
            }
        })
    });
    group.finish();
}

fn benchmark_geodtest_geodesic_inverse12(c: &mut Criterion) {
    let mut group = c.benchmark_group("Geodesic.inverse GeodTest adjusted");
    group.measurement_time(Duration::from_secs_f64(75.0));
    let data_arc = Arc::new(Mutex::new((
        util::as_vecs_num_sparse(&util::read_geodtest(), 0, &[0, 1, 3, 4]),
        geographiclib_rs::Geodesic::wgs84(),
    )));
    group.bench_function("Geodesic.inverse (rust impl, GeodTest.dat inputs)", |b| {
        let pair = data_arc.lock().unwrap();
        let data = &pair.0;
        let g = &pair.1;
        b.iter(|| {
            for line in data {
                let (_s12, _azi1, _azi2, _m12, _M12, _M21, _S12, _a12) =
                    g.inverse(line[0], line[1], line[2], line[3]);
            }
        })
    });
    group.finish();
}

// Benchmarks using logged cpp values as inputs: *_loggedcpp_*

// Note that using Arc and Mutex to access a shared dataset results
// in much lower overhead than cloning the dataset in each iteration.

// placeholder: Geodesic_A1m1f
// placeholder: Geodesic_A2m1f
// placeholder: Geodesic_A3coeff
// placeholder: Geodesic_A3f

// Note: For Geodesic_Astroid, see geomath_benchmark.rs

// placeholder: Geodesic_C1f
// placeholder: Geodesic_C1pf
// placeholder: Geodesic_C2f
// placeholder: Geodesic_C3coeff
// placeholder: Geodesic_C3f
// placeholder: Geodesic_C4coeff
// placeholder: Geodesic_C4f

fn benchmark_loggedcpp_geodesic_gen_direct(c: &mut Criterion) {
    // Format: this-in[_a _f] lat1 lon1 azi1 arcmode s12_a12 outmask result=a12 lat2-out lon2-out azi2-out s12-out m12-out M12-out M21-out S12-out
    let data = util::as_vecs_basic("Geodesic_GenDirect", 8);
    let geods = util::get_geodesic_lookup(&data);
    let data_arc = Arc::new(Mutex::new((data, geods)));
    c.bench_function("Geodesic._gen_direct (rust impl, cpp log inputs)", |b| {
        let pair = data_arc.lock().unwrap();
        let data = &pair.0;
        let geods = &pair.1;
        b.iter(|| {
            for line in data {
                let g = util::get_geodesic(line[0], line[1], &geods);
                let _result = g._gen_direct(line[2], line[3], line[4], line[5] != 0.0, line[6], line[7] as u64);
            }
        })
    });
}

// placeholder: Geodesic_GenDirectLine

fn benchmark_loggedcpp_geodesic_gen_inverse_azi(c: &mut Criterion) {
    // Format: this-in[_a _f] lat1 lon1 lat2 lon2 outmask result=a12 s12-out azi1-out azi2-out m12-out M12-out M21-out S12-out
    let data = util::as_vecs_basic("Geodesic_GenInverse_out7", 7);
    let geods = util::get_geodesic_lookup(&data);
    let data_arc = Arc::new(Mutex::new((data, geods)));
    c.bench_function("Geodesic._gen_inverse_azi (rust impl, cpp log inputs)", |b| {
        let pair = data_arc.lock().unwrap();
        let data = &pair.0;
        let geods = &pair.1;
        b.iter(|| {
            for line in data {
                let g = util::get_geodesic(line[0], line[1], &geods);
                let _result = g._gen_inverse_azi(line[2], line[3], line[4], line[5], line[6] as u64);
            }
        })
    });
}

fn benchmark_loggedcpp_geodesic_gen_inverse(c: &mut Criterion) {
    // Format: this-in[_a _f] lat1 lon1 lat2 lon2 outmask result=a12 s12-out salp1-out calp1-out salp2-out calp2-out m12-out M12-out M21-out S12-out
    let data = util::as_vecs_basic("Geodesic_GenInverse_out9", 7);
    let geods = util::get_geodesic_lookup(&data);
    let data_arc = Arc::new(Mutex::new((data, geods)));
    c.bench_function("Geodesic._gen_inverse (rust impl, cpp log inputs)", |b| {
        let pair = data_arc.lock().unwrap();
        let data = &pair.0;
        let geods = &pair.1;
        b.iter(|| {
            for line in data {
                let g = util::get_geodesic(line[0], line[1], &geods);
                let _result = g._gen_inverse(line[2], line[3], line[4], line[5], line[6] as u64);
            }
        })
    });
}

fn benchmark_loggedcpp_geodesic_new(c: &mut Criterion) {
    // Format: this-in[_a _f] lat1 lon1 lat2 lon2 outmask result=a12 s12-out salp1-out calp1-out salp2-out calp2-out m12-out M12-out M21-out S12-out
    let data_arc = Arc::new(Mutex::new(util::as_vecs_basic("Geodesic_Geodesic", 2)));
    c.bench_function("Geodesic.new (rust impl, cpp log inputs)", |b| {
        let data = data_arc.lock().unwrap();
        b.iter(|| {
            for i in 0..data.len() {
                let line = &data[i];
                let _result = geographiclib_rs::Geodesic::new(line[0], line[1]);
            }
        })
    });
}

// placeholder: Geodesic_InverseLine

fn benchmark_loggedcpp_geodesic_inverse_start(c: &mut Criterion) {
    // Format: this-in[_a _f] sbet1 cbet1 dn1 sbet2 cbet2 dn2 lam12 slam12 clam12 result=sig12 salp1-out calp1-out salp2-out calp2-out dnm-out
    let mut group = c.benchmark_group("Geodesic._InverseStart adjusted");
    group.measurement_time(Duration::from_secs_f64(7.5));

    let data = util::as_vecs_basic("Geodesic_InverseStart", 11);
    let geods = util::get_geodesic_lookup(&data);
    let data_arc = Arc::new(Mutex::new((data, geods)));
    group.bench_function("Geodesic._InverseStart (rust impl, cpp log inputs)", |b| {
        let pair = data_arc.lock().unwrap();
        let data = &pair.0;
        let geods = &pair.1;
        #[allow(non_snake_case)]
        let mut C1a: [f64; nC_] = [0.0 ; nC_];
        #[allow(non_snake_case)]
        let mut C2a: [f64; nC_] = [0.0 ; nC_];
        b.iter(|| {
            for line in data {
                let g = util::get_geodesic(line[0], line[1], &geods);
                let _result = g._InverseStart(line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], &mut C1a, &mut C2a);
            }
        })
    });

    group.finish();
}

fn benchmark_loggedcpp_geodesic_lambda12(c: &mut Criterion) {
    // Format: this-in[_a _f] sbet1 cbet1 dn1 sbet2 cbet2 dn2 salp1 calp1 slam120 clam120 diffp result=lam12 salp2-out calp2-out sig12-out ssig1-out csig1-out ssig2-out csig2-out eps-out domg12-out dlam12-out
    let data = util::as_vecs_basic("Geodesic_Lambda12", 13);
    let geods = util::get_geodesic_lookup(&data);
    let data_arc = Arc::new(Mutex::new((data, geods)));
    c.bench_function("Geodesic._Lambda12 (rust impl, cpp log inputs)", |b| {
        let pair = data_arc.lock().unwrap();
        let data = &pair.0;
        let geods = &pair.1;
        #[allow(non_snake_case)]
        let mut C1a: [f64; nC_] = [0.0 ; nC_];
        #[allow(non_snake_case)]
        let mut C2a: [f64; nC_] = [0.0 ; nC_];
        #[allow(non_snake_case)]
        let mut C3a: [f64; nC_] = [0.0 ; nC_];
        b.iter(|| {
            for line in data {
                let g = util::get_geodesic(line[0], line[1], &geods);
                let mut calp1 = line[9];
                let _result = g._Lambda12(line[2], line[3], line[4], line[5], line[6], line[7], line[8], &mut calp1, line[10], line[11], line[12] != 0.0, &mut C1a, &mut C2a, &mut C3a);
            }
        })
    });
}

fn benchmark_loggedcpp_geodesic_lengths(c: &mut Criterion) {
    // Format: this-in[_a _f] eps sig12 ssig1 csig1 dn1 ssig2 csig2 dn2 cbet1 cbet2 outmask s12b-out m12b-out m0-out M12-out M21-out
    let data = util::as_vecs_basic("Geodesic_Lengths", 13);
    let geods = util::get_geodesic_lookup(&data);
    let data_arc = Arc::new(Mutex::new((data, geods)));
    c.bench_function("Geodesic._Lengths (rust impl, cpp log inputs)", |b| {
        let pair = data_arc.lock().unwrap();
        let data = &pair.0;
        let geods = &pair.1;
        #[allow(non_snake_case)]
        let mut C1a: [f64; nC_] = [0.0 ; nC_];
        #[allow(non_snake_case)]
        let mut C2a: [f64; nC_] = [0.0 ; nC_];
        b.iter(|| {
            for line in data {
                let g = util::get_geodesic(line[0], line[1], &geods);
                let _result = g._Lengths(line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12] as u64, &mut C1a, &mut C2a);
            }
        })
    });
}

// Note: For Geodesic_SinCosSeries see geomath_benchmark.rs


criterion_group!(
    benches,
    geodesic_direct_benchmark,
    geodesic_inverse_benchmark,
    benchmark_geodtest_geodesic_direct12,
    benchmark_geodtest_geodesic_inverse12,
    benchmark_loggedcpp_geodesic_gen_direct,
    benchmark_loggedcpp_geodesic_gen_inverse_azi,
    benchmark_loggedcpp_geodesic_gen_inverse,
    benchmark_loggedcpp_geodesic_new,
    benchmark_loggedcpp_geodesic_inverse_start,
    benchmark_loggedcpp_geodesic_lambda12,
    benchmark_loggedcpp_geodesic_lengths,
);
criterion_main!(benches);
