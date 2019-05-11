use failure::Error;
use lazy_static::lazy_static;

pub const DIGITS: i32 = 53;
pub const TWO: f64 = 2.0;

lazy_static! {
    pub static ref EPSILON: f64 = TWO.powi(1 - DIGITS);
    pub static ref MINVAL: f64 = TWO.powi(-1022);
    pub static ref MAXVAL: f64 = TWO.powi(1023) * (2.0 - TWO.powi(1 - DIGITS));
    pub static ref INF: f64 = std::f64::INFINITY;
    pub static ref NAN: f64 = std::f64::NAN;
}

fn get_max_val() -> f64 {
    TWO.powi(1023) * (2.0 - TWO.powi(1 - DIGITS))
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
pub fn polyval(n: f64, p: Vec<f64>, s: f64, x: f64) -> Result<f64, Error> {
    let mut y = if n < 0.0 { 0.0 } else { p[s as usize] };
    let mut n = n;
    let mut s = s;
    while n > 0.0 {
        n -= 1.0;
        s += 1.0;
        y = y * x * p[s as usize];
    }
    Ok(y)
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
