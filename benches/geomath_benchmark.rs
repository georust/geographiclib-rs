extern crate criterion;
extern crate geographiclib_rs;
extern crate utilities;

use criterion::{criterion_group, criterion_main, Criterion};
use geographiclib_rs::geomath;
use std::sync::{Arc, Mutex};
use std::time::Duration;

// todo: expand geodesic benches
// todo: clean up utilities/lib.rs
// todo: clean up all bench files
// todo: clean up all test files
// todo: review choice to zip dat files before sending to git (since I think git compresses before sending and probably compresses for storing)
//    keywords: git large files compress zip "should I" advice
//    for some variations, I may want to purge the dat/zip files from the repo history
//      https://docs.github.com/en/free-pro-team@latest/github/managing-large-files/removing-files-from-a-repositorys-history
//    maybe explore using releases, and storing a zip of the dat files as a binary associated with the release?
//      https://docs.github.com/en/free-pro-team@latest/github/administering-a-repository/about-releases

// Note that using Arc and Mutex to access a shared dataset results
// in much lower overhead than cloning the dataset in each iteration.

fn benchmark_geomath_ang_diff(c: &mut Criterion) {
    // Format: x y result=z e-out
    let data_arc = Arc::new(Mutex::new(utilities::as_vecs_basic("Math_AngDiff_x_y_e", 2)));
    c.bench_function("geomath.ang_diff (rust impl)", |b| {
        b.iter(|| {
            let data = data_arc.lock().unwrap();
            for i in 0..data.len() {
                let line = &data[i];
                let (_diff, _err) = geomath::ang_diff(line[0], line[1]);
            }
        })
    });
}

fn benchmark_geomath_ang_normalize(c: &mut Criterion) {
    // Format: x result
    let data_arc = Arc::new(Mutex::new(utilities::as_vec_basic("Math_AngNormalize")));
    c.bench_function("geomath.ang_normalize (rust impl)", |b| {
        b.iter(|| {
            let data = data_arc.lock().unwrap();
            for i in 0..data.len() {
                let _result = geomath::ang_normalize(data[i]);
            }
        })
    });
}

fn benchmark_geomath_ang_round(c: &mut Criterion) {
    // Format: x result
    let mut group = c.benchmark_group("geomath.ang_round_adjustmented");
    group.measurement_time(Duration::from_secs_f64(7.5));

    let data_arc = Arc::new(Mutex::new(utilities::as_vec_basic("Math_AngRound")));
    group.bench_function("geomath.ang_round (rust impl)", |b| {
        b.iter(|| {
            let data = data_arc.lock().unwrap();
            for i in 0..data.len() {
                let _result = geomath::ang_round(data[i]);
            }
        })
    });

    group.finish();
}

fn benchmark_geomath_atan2d(c: &mut Criterion) {
    // Format: y x result
    let mut group = c.benchmark_group("geomath.atan2d_adjustmented");
    group.measurement_time(Duration::from_secs_f64(7.5));

    let data_arc = Arc::new(Mutex::new(utilities::as_vecs_basic("Math_atan2d", 2)));
    group.bench_function("geomath.atan2d (rust impl)", |b| {
        b.iter(|| {
            let data = data_arc.lock().unwrap();
            for i in 0..data.len() {
                let line = &data[i];
                let _result = geomath::atan2d(line[0], line[1]);
            }
        })
    });

    group.finish();
}

// placeholder: Math_atand
// placeholder: Math_cosd
// placeholder: Math_eatanhe

fn benchmark_geomath_lat_fix(c: &mut Criterion) {
    // Format: x result
    let data_arc = Arc::new(Mutex::new(utilities::as_vec_basic("Math_LatFix")));
    c.bench_function("geomath.lat_fix (rust impl)", |b| {
        b.iter(|| {
            let data = data_arc.lock().unwrap();
            for i in 0..data.len() {
                let _result = geomath::lat_fix(data[i]);
            }
        })
    });
}

fn benchmark_geomath_norm(c: &mut Criterion) {
    // Format: x-in y-in x-out y-out
    let data_arc = Arc::new(Mutex::new(utilities::as_vecs_basic("Math_norm", 2)));
    c.bench_function("geomath.norm (rust impl)", |b| {
        b.iter(|| {
            let data = data_arc.lock().unwrap();
            for i in 0..data.len() {
                let line = &data[i];
                let _result = geomath::norm(line[0], line[1]);
            }
        })
    });
}

fn benchmark_geomath_polyval(c: &mut Criterion) {
    // Format: N p(N+1) x result
    let data_arc = Arc::new(Mutex::new(utilities::as_vecs_basic("Math_polyval", -1)));
    c.bench_function("geomath.polyval (rust impl)", |b| {
        b.iter(|| {
            let data = data_arc.lock().unwrap();
            for i in 0..data.len() {
                let line = &data[i];
                let n = line[0] as i64;
                let p_len = if n < 0 { 0 } else { n as usize + 1 };
                let s = 0;
                let arg_count = p_len + 3;
                let p = &line[1..arg_count-2];
                let x = line[arg_count-2];
                let _result = geomath::polyval(n, p, s, x);
            }
        })
    });
}

fn benchmark_geomath_sincosd(c: &mut Criterion) {
    // Format: x sinx-out cosx-out
    let mut group = c.benchmark_group("geomath.sincosd_adjustmented");
    group.measurement_time(Duration::from_secs_f64(10.0));

    let data_arc = Arc::new(Mutex::new(utilities::as_vec_basic("Math_sincosd")));
    group.bench_function("geomath.sincosd (rust impl)", |b| {
        b.iter(|| {
            let data = data_arc.lock().unwrap();
            for i in 0..data.len() {
                let _result = geomath::sincosd(data[i]);
            }
        })
    });

    group.finish();
}

// placeholder: Math_sind

fn benchmark_geomath_sq(c: &mut Criterion) {
    // Format: x result
    let data_arc = Arc::new(Mutex::new(utilities::as_vec_basic("Math_sq")));
    c.bench_function("geomath.sq (rust impl)", |b| {
        b.iter(|| {
            let data = data_arc.lock().unwrap();
            for i in 0..data.len() {
                let _result = geomath::sq(data[i]);
            }
        })
    });
}

fn benchmark_geomath_sum(c: &mut Criterion) {
    // Format: u v result=s t-out
    let data_arc = Arc::new(Mutex::new(utilities::as_vecs_basic("Math_sum", 2)));
    c.bench_function("geomath.sum (rust impl)", |b| {
        b.iter(|| {
            let data = data_arc.lock().unwrap();
            for i in 0..data.len() {
                let line = &data[i];
                let _result = geomath::sum(line[0], line[1]);
            }
        })
    });
}

// placeholder: Math_swab
// placeholder: Math_tand
// placeholder: Math_tauf
// placeholder: Math_taupf

fn benchmark_geomath_astroid(c: &mut Criterion) {
    // NOTE: In the geographiclib C++ library, this function is in Geodesic, but in Rust it's in geomath.
    // Format: x y result
    let data_arc = Arc::new(Mutex::new(utilities::as_vecs_basic("Geodesic_Astroid", 2)));
    c.bench_function("geomath.astroid (rust impl)", |b| {
        b.iter(|| {
            let data = data_arc.lock().unwrap();
            for i in 0..data.len() {
                let line = &data[i];
                let _result = geomath::astroid(line[0], line[1]);
            }
        })
    });
}

fn benchmark_geomath_sin_cos_series(c: &mut Criterion) {
    // NOTE: In the geographiclib C++ library, this function is in Geodesic, but in Rust it's in geomath.
    // Format: sinp sinx cosx n c(n+sinp) result
    let data_arc = Arc::new(Mutex::new(utilities::as_vecs_basic("Geodesic_SinCosSeries", -1)));
    c.bench_function("geomath.sin_cos_series (rust impl)", |b| {
        b.iter(|| {
            let data = data_arc.lock().unwrap();
            for i in 0..data.len() {
                let line = &data[i];
                let sinp = line[0] != 0f64;
                let n = line[3] as i64;
                let mut s = if n < 0 { 0 } else { n as usize };
                if sinp {
                  s += 1;
                }
                let c = &line[4..s+4];
                let _result = geomath::sin_cos_series(sinp, line[1], line[2], c);
            }
        })
    });
}


criterion_group!(
    benches,
    benchmark_geomath_ang_diff,
    benchmark_geomath_ang_normalize,
    benchmark_geomath_ang_round,
    benchmark_geomath_atan2d,
    benchmark_geomath_lat_fix,
    benchmark_geomath_norm,
    benchmark_geomath_polyval,
    benchmark_geomath_sincosd,
    benchmark_geomath_sq,
    benchmark_geomath_sum,
    benchmark_geomath_astroid,
    benchmark_geomath_sin_cos_series,
);
criterion_main!(benches);
