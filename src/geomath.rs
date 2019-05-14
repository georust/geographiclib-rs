use lazy_static::lazy_static;

pub const DIGITS: u64 = 53;
pub const TWO: f64 = 2.0;

lazy_static! {
    pub static ref EPSILON: f64 = TWO.powi(1 - DIGITS as i32);
    pub static ref MINVAL: f64 = TWO.powi(-1022);
    pub static ref MAXVAL: f64 = TWO.powi(1023) * (2.0 - TWO.powi(1 - DIGITS as i32));
    pub static ref INF: f64 = std::f64::INFINITY;
    pub static ref NAN: f64 = std::f64::NAN;
}

pub fn get_epsilon() -> f64 {
    TWO.powi(1 - DIGITS as i32)
}

pub fn get_min_val() -> f64 {
    TWO.powi(1023) * (2.0 - TWO.powi(1 - DIGITS as i32))
}

pub fn get_max_val() -> f64 {
    TWO.powi(1023) * (2.0 - TWO.powi(1 - DIGITS as i32))
}

// Square
pub fn sq(x: f64) -> f64 {
    x.powi(2)
}

// Real cube root
pub fn cbrt(x: f64) -> f64 {
    let y = x.abs().powf(1.0 / 3.0);
    if x >= 0.0 {
        y
    } else {
        -y
    }
}

// Normalize a two-vector
pub fn norm(x: f64, y: f64) -> (f64, f64) {
    let r = x.hypot(y);
    (x / r, y / r)
}

// Error free transformation of a sum
pub fn sum(u: f64, v: f64) -> (f64, f64) {
    let s = u + v;
    let up = s - v;
    let vpp = s - up;
    let up = up - u;
    let vpp = vpp - v;
    let t = -(up + vpp);
    (s, t)
}

// Evaluate a polynomial
pub fn polyval(n: i64, p: &Vec<f64>, s: usize, x: f64) -> f64 {
    let mut s = s;
    let mut n = n;
    let mut y = if n < 0 { 0.0 } else { p[s] };
    while n > 0 {
        n -= 1;
        s = s.checked_add(1).expect("");
        y = y * x + p[s];
    }
    y
}

// Round an angle so taht small values underflow to 0
pub fn ang_round(x: f64) -> f64 {
    // The makes the smallest gap in x = 1/16 - nextafter(1/16, 0) = 1/2^57
    // for reals = 0.7 pm on the earth if x is an angle in degrees.  (This
    // is about 1000 times more resolution than we get with angles around 90
    // degrees.)  We use this to avoid having to deal with near singular
    // cases when x is non-zero but tiny (e.g., 1.0e-200).
    let z = 1.0 / 16.0;
    let mut y = x.abs();
    // The compiler mustn't "simplify" z - (z - y) to y
    if y < z {
        y = z - (z - y);
    };
    if x == 0.0 {
        0.0
    } else {
        if x < 0.0 {
            -y
        } else {
            y
        }
    }
}

// reduce angle to (-180,180]
pub fn ang_normalize(x: f64) -> f64 {
    let mut y = x % 2.0 * std::f64::consts::PI;
    if x == 0.0 {
        y = x;
    };
    if y <= -180.0 {
        y + 360.0
    } else {
        if y <= 180.0 {
            y
        } else {
            y - 360.0
        }
    }
}

// Replace angles outside [-90,90] with NaN
pub fn lat_fix(x: f64) -> f64 {
    if x.abs() > 90.0 {
        std::f64::NAN
    } else {
        x
    }
}

// compute y - x and reduce to [-180,180] accurately
pub fn ang_diff(x: f64, y: f64) -> (f64, f64) {
    let (d, t) = sum(ang_normalize(-x), ang_normalize(y));
    let d = ang_normalize(d);
    if d == 180.0 && t > 0.0 {
        sum(-180.0, t)
    } else {
        sum(d, t)
    }
}

// Compute sine and cosine of x in degrees
pub fn sincosd(x: f64) -> (f64, f64) {
    let mut r = x % 2.0 * std::f64::consts::PI;
    let mut q = if r.is_nan() {
        std::f64::NAN
    } else {
        (r / 90.0 + 0.5).floor()
    };
    r -= 90.0 * q;
    r = r.to_radians();
    let mut s = r.sin();
    let mut c = r.cos();
    q = q % 4.0;
    if q == 1.0 {
        let _s = s;
        s = c;
        c = -_s;
    } else if q == 2.0 {
        s = -s;
        c = -c;
    } else if q == 3.0 {
        let _s = s;
        s = -c;
        c = _s
    }
    if x == 0.0 {
        s = x;
    } else {
        s = 0.0 + s;
        c = 0.0 + c;
    }
    (s, c)
}

// Compute atan2(y, x) with result in degrees
pub fn atan2d(y: f64, x: f64) -> f64 {
    let mut x = x;
    let mut y = y;
    let mut q = if y.abs() > x.abs() {
        let _x = x;
        x = y;
        y = _x;
        2.0
    } else {
        0.0
    };
    if x < 0.0 {
        q += 1.0;
        x = -x;
    }
    let mut ang = y.atan2(x).to_degrees();
    if q == 1.0 {
        ang = if y >= 0.0 { 180.0 - ang } else { -180.0 - ang };
    } else if q == 2.0 {
        ang = 90.0 - ang;
    } else if q == 3.0 {
        ang = -90.0 + ang;
    }
    ang
}

// test for finitness
pub fn isfinite(x: f64) -> bool {
    x.abs() <= get_max_val()
}

// Functions that used to be inside Geodesic
pub fn sin_cos_series(sinp: bool, sinx: f64, cosx: f64, c: Vec<i64>) -> f64 {
    let mut k = c.len();
    let mut n = k - if sinp { 1 } else { 0 };
    let ar: f64 = 2.0 * (cosx - sinx) * (cosx + sinx);
    let mut y1 = 0.0;
    let mut y0: f64 = if n != 0 {
        k -= 1;
        c[k] as f64
    } else {
        0.0
    };
    n = n / 2;
    while n > 0 {
        n -= 1;
        k -= 1;
        y1 = ar * y0 - y1 + c[k] as f64;
        k -= 1;
        y0 = ar * y1 - y0 + c[k] as f64;
    }
    2.0 * sinx * cosx * if sinp { y0 } else { cosx * (y0 - y1) }
}

// Solve stroid equation
pub fn astroid(x: f64, y: f64) -> f64 {
    let p = sq(x);
    let q = sq(y);
    let r = (p + q - 1.0) / 6.0;
    if !(q == 0.0 && r <= 0.0) {
        let s = p * q / 4.0;
        let r2 = sq(r);
        let r3 = r * r2;
        let disc = s * (s + 2.0 * r3);
        let mut u = r;
        if disc >= 0.0 {
            let mut t3 = s + r3;
            t3 += if t3 < 0.0 { -disc.sqrt() } else { disc.sqrt() };
            let t = t3.cbrt();
            u += t + if t != 0.0 { r2 / t } else { 0.0 };
        } else {
            let ang = (-disc).sqrt().atan2(-(s + r3));
            u += 2.0 * r * (ang / 3.0).cos();
        }
        let v = (sq(u) + q).sqrt();
        let uv = if u < 0.0 { q / (v - u) } else { u + v };
        let w = (uv - q) / (2.0 * v);
        uv / ((uv + sq(w)).sqrt() + w)
    } else {
        0.0
    }
}

pub fn _A1m1f(eps: f64, geodesic_order: i64) -> f64 {
    let coeff = vec![1.0, 4.0, 64.0, 0.0, 256.0];
    let m: i64 = geodesic_order / 2;
    let t = polyval(m, &coeff, 0, sq(eps)) / coeff[(m + 1) as usize] as f64;
    (t + eps) / (1.0 - eps)
}

pub fn _C1f(eps: f64, c: &mut Vec<f64>, geodesic_order: i64) {
    let coeff = vec![
        -1.0, 6.0, -16.0, 32.0, -9.0, 64.0, -128.0, 2048.0, 9.0, -16.0, 768.0, 3.0, -5.0, 512.0,
        -7.0, 1280.0, -7.0, 2048.0,
    ];
    let eps2 = sq(eps);
    let mut d = eps;
    let mut o = 0;
    for l in 1..=geodesic_order {
        let m = ((geodesic_order - l) / 2) as i64;
        c[l as usize] =
            d * polyval(m, &coeff, o as usize, eps2) / coeff[(o + m + 1) as usize] as f64;
        o += m + 2;
        d *= eps;
    }
}

pub fn _C1pf(eps: f64, c: &mut Vec<f64>, geodesic_order: i64) {
    let coeff = vec![
        205.0, -432.0, 768.0, 1536.0, 4005.0, -4736.0, 3840.0, 12288.0, -225.0, 116.0, 384.0,
        -7173.0, 2695.0, 7680.0, 3467.0, 7680.0, 38081.0, 61440.0,
    ];
    let eps2 = sq(eps);
    let mut d = eps;
    let mut o = 0;
    for l in 1..=geodesic_order {
        let m = (geodesic_order - l) / 2;
        c[l as usize] =
            d * polyval(m as i64, &coeff, o as usize, eps2) / coeff[(o + m + 1) as usize] as f64;
        o += m + 2;
        d *= eps;
    }
}

pub fn _A2m1f(eps: f64, geodesic_order: i64) -> f64 {
    let coeff = vec![-11.0, -28.0, -192.0, 0.0, 256.0];
    let m: i64 = geodesic_order / 2;
    let t = polyval(m, &coeff, 0, sq(eps)) / coeff[(m + 1) as usize] as f64;
    (t - eps) / (1.0 + eps)
}

pub fn _C2f(eps: f64, c: &mut Vec<f64>, geodesic_order: i64) {
    let coeff = vec![
        1.0, 2.0, 16.0, 32.0, 35.0, 64.0, 384.0, 2048.0, 15.0, 80.0, 768.0, 7.0, 35.0, 512.0, 63.0,
        1280.0, 77.0, 2048.0,
    ];
    let eps2 = sq(eps);
    let mut d = eps;
    let mut o = 0;
    for l in 1..=geodesic_order {
        let m = (geodesic_order - l) / 2;
        c[l as usize] =
            d * polyval(m as i64, &coeff, o as usize, eps2) / coeff[(o + m + 1) as usize] as f64;
        o += m + 2;
        d *= eps;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    // Results for the assertions are taken by running the python implementation

    #[test]
    fn test__C2f() {
        let mut c = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        _C2f(0.12, &mut c, 6);
        assert_eq!(
            c,
            vec![
                1.0,
                0.0601087776,
                0.00270653103,
                0.000180486,
                1.4215824e-05,
                1.22472e-06,
                1.12266e-07
            ]
        )
    }

    #[test]
    fn test__A2m1f() {
        assert_eq!(_A2m1f(0.12, 6), -0.11680607884285714);
    }

    #[test]
    fn test__C1pf() {
        let mut c = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        _C1pf(0.12, &mut c, 6);
        assert_eq!(
            c,
            vec![
                1.0,
                0.059517321000000005,
                0.004421053215,
                0.0005074200000000001,
                6.997613759999999e-05,
                1.1233080000000001e-05,
                1.8507366e-06
            ]
        )
    }

    #[test]
    fn test__C1f() {
        let mut c = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        _C1f(0.12, &mut c, 6);
        assert_eq!(
            c,
            vec![
                1.0,
                -0.059676777599999994,
                -0.000893533122,
                -3.57084e-05,
                -2.007504e-06,
                -1.3607999999999999e-07,
                -1.0205999999999999e-08
            ]
        )
    }

    #[test]
    fn test__A1m1f() {
        assert_eq!(_A1m1f(0.12, 6), 0.1404582405272727);
    }

    #[test]
    fn test_astroid() {
        assert_eq!(astroid(21.0, 12.0), 23.44475767500982);
    }

    #[test]
    fn test_sin_cos_series() {
        assert_eq!(
            sin_cos_series(true, 0.12, 0.21, vec![1, 2]),
            0.10079999999999999
        );
    }
}
