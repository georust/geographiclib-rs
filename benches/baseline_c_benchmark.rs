extern crate criterion;
extern crate geographiclib;

mod util;

use criterion::{criterion_group, criterion_main, Criterion};
use std::time::Duration;

fn geodesic_direct_baseline(c: &mut Criterion) {
    let (mode, records) = util::load_test_data();

    let mut group = c.benchmark_group("direct (c wrapper)");
    if mode == util::TEST_MODE_FULL {
        group.measurement_time(Duration::from_secs(30));
    }
    group.bench_function(mode, |b| {
        let geod = geographiclib::Geodesic::wgs84();
        b.iter(|| {
            for r in &records {
                let (_lat2, _lon2, _azi2) = geod.direct(r.lat1, r.lon1, r.azi1, r.s12);
            }
        })
    });
    group.finish();
}

fn geodesic_inverse_baseline(c: &mut Criterion) {
    let (mode, records) = util::load_test_data();

    let mut group = c.benchmark_group("inverse (c wrapper)");
    if mode == util::TEST_MODE_FULL {
        group.measurement_time(Duration::from_secs(50));
    }
    group.bench_function(mode, |b| {
        let geod = geographiclib::Geodesic::wgs84();
        b.iter(|| {
            for r in &records {
                let (_s12, _azi1, _azi2, _a12) = geod.inverse(r.lat1, r.lon1, r.lat2, r.lon2);
            }
        })
    });
    group.finish();
}

criterion_group!(benches, geodesic_direct_baseline, geodesic_inverse_baseline,);
criterion_main!(benches);
