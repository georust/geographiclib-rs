// Integration tests related to geomath.

extern crate utilities;

use geographiclib_rs::geomath;
use utilities::{assert_delta, test_basic};

// reminder: use "cargo test -- --nocapture" to let tests print to console.
// reminder: if a test fails especially badly, review the C++ "instrumented-crude" logging for possible bugs there

// todo: decide between keeping these as integration tests vs unit tests vs mixed...
//       why integration: external files and parsing code (not self-contained)
//       why unit: sometimes want access to private/internal data (e.g. consts, constructors)

#[test]
fn test_ang_diff_vs_cpp() {
    // Format: x y result=z e-out
    test_basic("Math_AngDiff_x_y_e", 4, |line_num, items| {
        let result = geomath::ang_diff(items[0], items[1]);
        // todo: try to pass allow_sign_change=false
        assert_delta!(items[2], result.0, 0.0, true, "result (result.0)", line_num);
        // todo: try to pass allow_sign_change=false
        assert_delta!(items[3], result.1, 0.0, true, "e (result.1)", line_num);
    });
}

#[test]
fn test_ang_normalize_vs_cpp() {
    // Format: x result
    test_basic("Math_AngNormalize", 2, |line_num, items| {
        let result = geomath::ang_normalize(items[0]);
        // todo: try to pass allow_sign_change=false
        assert_delta!(items[1], result, 0.0, true, "result", line_num);
    });
}

#[test]
fn test_ang_round_vs_cpp() {
    // Format: x result
    test_basic("Math_AngRound", 2, |line_num, items| {
        let result = geomath::ang_round(items[0]);
        assert_delta!(items[1], result, 0.0, false, "result", line_num);
    });
}

#[test]
fn test_atan2d_vs_cpp() {
    // Format: y x result
    test_basic("Math_atan2d", 3, |line_num, items| {
        let result = geomath::atan2d(items[0], items[1]);
        // todo: try to reduce allowed difference to 0.0
        assert_delta!(items[2], result, 5e-16, false, "result", line_num);
    });
}

// #[test]
// fn test_atand_vs_cpp() {
//     // Format: x result
//     test_basic("Math_atand", 2, |line_num, items| {
//         let result = geomath::atand(items[0]);
//         assert_delta!(items[1], result, 0.0, false, "result", line_num);
//     });
// }

#[test]
fn test_consts_vs_cpp() {
    // Format: digits digits10 extra_digits bigendian pi degree GEOGRAPHICLIB_PRECISION GEOGRAPHICLIB_WORDS_BIGENDIAN
    // todo: review commented-out tests to see if there's an equivalent const or simple calculation we can compare
    let items = utilities::read_consts_basic("Math_consts", 8);
    let line_num = 2; // todo: remove once assert_delta stops requiring a line number
    assert_delta!(items[0], geomath::DIGITS as f64, 0.0, false, "DIGITS", line_num);
    // assert_delta!(items[1], geomath::DIGITS10 as f64, 0.0, false, "DIGITS10", line_num);
    // assert_delta!(items[2], geomath::EXTRA_DIGITS, 0.0, false, "EXTRA_DIGITS", line_num);
    // assert_delta!(items[3], geomath::BIGENDIAN, 0.0, false, "BIGENDIAN", line_num);
    assert_delta!(items[4], std::f64::consts::PI, 0.0, false, "PI", line_num);
    // assert_delta!(items[5], geomath::DEGREE, 0.0, false, "DEGREE", line_num);
    // assert_delta!(items[6], geomath::GEOGRAPHICLIB_PRECISION, 0.0, false, "GEOGRAPHICLIB_PRECISION", line_num);
    // assert_delta!(items[7], geomath::GEOGRAPHICLIB_WORDS_BIGENDIAN, 0.0, false, "GEOGRAPHICLIB_WORDS_BIGENDIAN", line_num);
}

// #[test]
// fn test_cosd_vs_cpp() {
//     // Format: x result
//     test_basic("Math_cosd", 2, |line_num, items| {
//         let result = geomath::cosd(items[0]);
//         assert_delta!(items[1], result, 0.0, false, "result", line_num);
//     });
// }

// #[test]
// fn test_eatanhe_vs_cpp() {
//     // Format: x es result
//     test_basic("Math_eatanhe", 3, |line_num, items| {
//         let result = geomath::eatanhe(items[0], items[1]);
//         assert_delta!(items[2], result, 0.0, false, "result", line_num);
//     });
// }

#[test]
fn test_lat_fix_vs_cpp() {
    // Format: x result
    test_basic("Math_LatFix", 2, |line_num, items| {
        let result = geomath::lat_fix(items[0]);
        assert_delta!(items[1], result, 0.0, false, "result", line_num);
    });
}

#[test]
fn test_norm_vs_cpp() {
    // Format: x-in y-in x-out y-out
    test_basic("Math_norm", 4, |line_num, items| {
        let result = geomath::norm(items[0], items[1]);
        assert_delta!(items[2], result.0, 0.0, false, "x-out (result.0)", line_num);
        assert_delta!(items[3], result.1, 0.0, false, "y-out (result.1)", line_num);
    });
}

#[test]
fn test_polyval_vs_cpp() {
    // Format: N p(N+1) x result
    test_basic("Math_polyval", -1, |line_num, items| {
        assert!(items.len() > 2, "Expected a minimum of 3 items per line. Line {} had {}.", line_num, items.len());
        let n = items[0] as i64;
        let p_len = if n < 0 { 0 } else { n as usize + 1 };
        let arg_count = p_len + 3;
        assert!(items.len() == arg_count, "Expected {} items on line {}, based on first item, but found {}.", arg_count, line_num, items.len());
        let p = &items[1..arg_count-2];
        let x = items[arg_count-2];
        assert_eq!(p_len, p.len(), "Internal error: On line {}, tried to construct slice size {} but found {}", line_num, p_len, p.len());
        let result = geomath::polyval(n, p, 0, x);
        assert_delta!(items[arg_count - 1], result, 0.0, false, "result", line_num);
    });
}

#[test]
fn test_sincosd_vs_cpp() {
    // Format: x sinx-out cosx-out
    test_basic("Math_sincosd", 3, |line_num, items| {
        let result = geomath::sincosd(items[0]);
        // todo: try to reduce allowed difference to 0.0
        // todo: try to pass allow_sign_change=false
        assert_delta!(items[1], result.0, 2e-16, true, "sinx-out (result.0)", line_num);
        // todo: try to reduce allowed difference to 0.0
        // todo: try to pass allow_sign_change=false
        assert_delta!(items[2], result.1, 2e-16, true, "cosx-out (result.1)", line_num);
    });
}

// #[test]
// fn test_sind_vs_cpp() {
//     // Format: x result
//     test_basic("Math_sind", 2, |line_num, items| {
//         let result = geomath::sind(items[0]);
//         assert_delta!(items[1], result, 0.0, false, "result", line_num);
//     });
// }

#[test]
fn test_sq_vs_cpp() {
    // Format: x result
    test_basic("Math_sq", 2, |line_num, items| {
        let result = geomath::sq(items[0]);
        assert_delta!(items[1], result, 0.0, false, "result", line_num);
    });
}

#[test]
fn test_sum_vs_cpp() {
    // Format: u v result=s t-out
    test_basic("Math_sum", 4, |line_num, items| {
        let result = geomath::sum(items[0], items[1]);
        assert_delta!(items[2], result.0, 0.0, false, "result (result.0)", line_num);
        assert_delta!(items[3], result.1, 0.0, false, "t-out (result.1)", line_num);
    });
}

// #[test]
// fn test_swab_vs_cpp() {
//     // Format: x result
//     test_basic("Math_swab", 2, |line_num, items| {
//         let result = geomath::swab(items[0] as isize);
//         assert_delta!(items[1], result as f64, 0.0, false, "result", line_num);
//     });
// }

// #[test]
// fn test_tand_vs_cpp() {
//     // Format: x result
//     test_basic("Math_tand", 2, |line_num, items| {
//         let result = geomath::tand(items[0]);
//         assert_delta!(items[1], result, 0.0, false, "result", line_num);
//     });
// }

// #[test]
// fn test_tauf_vs_cpp() {
//     // Format: taup es result
//     test_basic("Math_tauf", 3, |line_num, items| {
//         let result = geomath::tauf(items[0], items[1]);
//         assert_delta!(items[2], result, 0.0, false, "result", line_num);
//     });
// }

// #[test]
// fn test_taupf_vs_cpp() {
//     // Format: tau es result
//     test_basic("Math_taupf", 3, |line_num, items| {
//         let result = geomath::tauof(items[0], items[1]);
//         assert_delta!(items[2], result, 0.0, false, "result", line_num);
//     });
// }

#[test]
fn test_astroid_vs_cpp() {
    // Format: x y result
    // Note: In the geographiclib C++ library, this function is in Geodesic, but in Rust it's in geomath.
    test_basic("Geodesic_Astroid", 3, |line_num, items| {
        let result = geomath::astroid(items[0], items[1]);
        // todo: try to reduce allowed difference to 0.0
        // todo: try to pass allow_sign_change=false
        assert_delta!(items[2], result, 1e-15, true, "result", line_num);
    });
}

#[test]
fn test_sin_cos_series_vs_cpp() {
    // Format: sinp sinx cosx n c(n+sinp) result
    // Note: In the geographiclib C++ library, this function is in Geodesic, but in Rust it's in geomath.
    test_basic("Geodesic_SinCosSeries", -1, |line_num, items| {
        assert!(items.len() > 4, "Expected a minimum of 5 items per line. Line {} had {}.", line_num, items.len());
        let sinp = items[0] != 0f64;
        let n = items[3] as i64;
        let mut s = if n < 0 { 0 } else { n as usize };
        if sinp {
            s += 1;
        }
        let arg_count = s + 5;
        assert!(items.len() == arg_count, "Expected {} items on line {}, based on first item, but found {}.", arg_count, line_num, items.len());
        let c = &items[4..s+4];
        assert_eq!(s, c.len(), "Internal error: On line {}, tried to construct slice size {} but found {}", line_num, s, c.len());
        let result = geomath::sin_cos_series(sinp, items[1], items[2], c);
        assert_delta!(items[arg_count - 1], result, 0.0, false, "result", line_num);
    });
}
