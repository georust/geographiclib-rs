extern crate criterion;
extern crate geographiclib_rs;
extern crate utilities;

use criterion::{criterion_group, criterion_main, Criterion};
use geographiclib_rs::geomath;
use std::fs::File;

fn geomath_benchmark_ang_diff(c: &mut Criterion) {
    let file = File::open("test_fixtures/Math_AngDiff_x_y_e.dat").unwrap();
    let reqs = utilities::DataRequirements {
        op_type: "Math_AngDiff_x_y_e",
        input_sizes: &[3],
        output_sizes: &[2],
        non_nums_in: utilities::NON_NUMS_NONE,
        non_nums_out: utilities::NON_NUMS_NONE,
    };
    let inputs_num_all: Vec<Vec<f64>> = reqs.read_data_lines(&file).map(|parsed_line| parsed_line.inputs_num).collect();

    c.bench_function("ang_diff (rust impl)", |b| {
        b.iter(|| {
            for inputs_num in inputs_num_all.clone() {
                let (_diff, _err) = geomath::ang_diff(inputs_num[0], inputs_num[1]);
            }
        })
    });
}

criterion_group!(
    benches,
    geomath_benchmark_ang_diff,
);
criterion_main!(benches);
