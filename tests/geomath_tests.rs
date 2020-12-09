// Integration tests related to geomath.

extern crate utilities;

use std::fs::File;
use geographiclib_rs::geomath;
use utilities::assert_delta;

// reminder: use "cargo test -- --nocapture" to let tests print to console.

#[test]
fn test_ang_diff() {
  // Format: x y e ; ; e result ; testTimeSeconds
  let file = File::open("test_fixtures/Math_AngDiff_x_y_e.dat").unwrap();
  let reqs = utilities::DataRequirements {
    op_type: "Math_AngDiff_x_y_e",
    input_sizes: &[3],
    output_sizes: &[2],
    non_nums_in: utilities::NON_NUMS_NONE,
    non_nums_out: utilities::NON_NUMS_NONE,
  };
  let iter = reqs.read_data_lines(&file);
  for parsed_line in iter {
    let x = parsed_line.inputs_num[0];
    let y = parsed_line.inputs_num[1];
    let result = geomath::ang_diff(x, y);
    assert_delta!(parsed_line.outputs_num[1], result.0, 0.0, false, "result (result.0)");
    assert_delta!(parsed_line.outputs_num[0], result.1, 0.0, false, "e (result.1)");
  }
}