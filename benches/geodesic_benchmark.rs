extern crate criterion;
extern crate geographiclib_rs;

mod util;

use criterion::{criterion_group, criterion_main, Criterion};
use std::time::Duration;

use geographiclib_rs::{DirectGeodesic, InverseGeodesic, PolygonArea, Winding};

fn geodesic_direct_benchmark(c: &mut Criterion) {
    let (mode, records) = util::load_test_data();

    let mut group = c.benchmark_group("direct (rust impl)");
    if mode == util::TEST_MODE_FULL {
        group.measurement_time(Duration::from_secs(35));
    }
    group.bench_function(mode, |b| {
        let geod = geographiclib_rs::Geodesic::wgs84();
        b.iter(|| {
            for r in &records {
                let (_lat2, _lon2, _azi2) = geod.direct(r.lat1, r.lon1, r.azi1, r.s12);
            }
        })
    });
    group.finish();
}

fn geodesic_inverse_benchmark(c: &mut Criterion) {
    let (mode, records) = util::load_test_data();

    let mut group = c.benchmark_group("inverse (rust impl)");
    if mode == util::TEST_MODE_FULL {
        group.measurement_time(Duration::from_secs(70));
    }
    group.bench_function(mode, |b| {
        let geod = geographiclib_rs::Geodesic::wgs84();
        b.iter(|| {
            for r in &records {
                let (_s12, _azi1, _azi2, _a12) = geod.inverse(r.lat1, r.lon1, r.lat2, r.lon2);
            }
        })
    });
    group.finish();
}

fn polygon_area_benchmark(c: &mut Criterion) {
    let geojson_str = std::fs::read_to_string("test_fixtures/countries.geojson").unwrap();
    let geojson: geojson::GeoJson = geojson_str.parse().unwrap();
    let features = match geojson {
        geojson::GeoJson::FeatureCollection(fc) => fc.features,
        _ => panic!("expected FeatureCollection"),
    };

    // Extract all polygon rings as Vec<Vec<(lat, lon)>>
    let mut rings: Vec<Vec<(f64, f64)>> = Vec::new();
    for feature in &features {
        let geometry = feature.geometry.as_ref().unwrap();
        match &geometry.value {
            geojson::Value::Polygon(coords) => {
                for ring in coords {
                    rings.push(ring.iter().map(|c| (c[1], c[0])).collect());
                }
            }
            geojson::Value::MultiPolygon(polys) => {
                for polygon in polys {
                    for ring in polygon {
                        rings.push(ring.iter().map(|c| (c[1], c[0])).collect());
                    }
                }
            }
            _ => {}
        }
    }

    let mut group = c.benchmark_group("polygon_area");
    group.bench_function("countries", |b| {
        let geod = geographiclib_rs::Geodesic::wgs84();
        b.iter(|| {
            let mut area_sum = 0.0;
            for ring in &rings {
                let mut pa = PolygonArea::new(&geod, Winding::Clockwise);
                for &(lat, lon) in ring {
                    pa.add_point(lat, lon);
                }
                let (_perimeter, area, _point_count) = pa.compute(true);
                area_sum += area;
            }
            // Sanity check, some light googling confirms the land on earth is indeed about 150_000_000_000_000
            // so this is order of magnitude correct.
            assert_eq!(area_sum.round(), 147_362_836_335_330.0);
        })
    });
    group.finish();
}

criterion_group!(
    benches,
    geodesic_direct_benchmark,
    geodesic_inverse_benchmark,
    polygon_area_benchmark,
);
criterion_main!(benches);
