//! Geocentric calculations.
//!
//! Geocentric coordinates are also known as earth-centered, earth-fixed (ECEF)
//! coordinates. The origin is at the center of the earth. The Z axis passes
//! through the north pole, the X axis passes through latitude 0 and longitude
//! 0, and the Y axis completes a right-handed coordinate system.

use core::f64::consts::FRAC_PI_6;
use std::sync;

use crate::geodesic::{WGS84_A, WGS84_F};
use crate::geomath;

/// A converter between geodetic coordinates and geocentric coordinates.
///
/// This follows GeographicLib's C++ `Geocentric` implementation, which adapts
/// Vermeille's 2002 reverse transform with Karney's robustness improvements.
/// The general reverse branch uses Vermeille's 2011 formulas for the cubic
/// solution while preserving GeographicLib's special-case handling.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Geocentric {
    a: f64,
    f: f64,
    e2: f64,
    e2m: f64,
    e2a: f64,
    e4a: f64,
    maxrad: f64,
}

static WGS84_GEOCENTRIC: sync::OnceLock<Geocentric> = sync::OnceLock::new();

struct ReverseResult {
    lat: f64,
    lon: f64,
    height: f64,
    sphi: f64,
    cphi: f64,
    slam: f64,
    clam: f64,
}

impl Geocentric {
    /// Create a geocentric converter for the WGS 84 ellipsoid.
    #[inline]
    pub fn wgs84() -> Self {
        *WGS84_GEOCENTRIC.get_or_init(|| Self::new_unchecked(WGS84_A, WGS84_F))
    }

    /// Create a geocentric converter for an ellipsoid of revolution.
    ///
    /// Setting `f = 0` gives a sphere; negative `f` gives a prolate ellipsoid.
    ///
    /// # Arguments
    ///
    /// - `a` - equatorial radius (meters).
    /// - `f` - flattening of the ellipsoid.
    ///
    /// # Panics
    /// Panics if `a` is not positive and finite, or if `(1 - f) * a` is not
    /// positive and finite.
    #[inline]
    pub fn new(a: f64, f: f64) -> Self {
        assert!(a.is_finite() && a > 0.0, "semi-major axis must be positive");
        assert!(
            f.is_finite() && ((1. - f) * a).is_finite() && (1. - f) * a > 0.,
            "polar semi-axis must be positive"
        );
        Self::new_unchecked(a, f)
    }

    #[inline]
    const fn new_unchecked(a: f64, f: f64) -> Self {
        let e2 = f * (2.0 - f);
        let e2m = (1.0 - f) * (1.0 - f);
        Self {
            a,
            f,
            e2,
            e2m,
            // Inline of `e2.abs()`; `f64::abs` is not yet usable from a `const fn`.
            e2a: if e2 < 0.0 { -e2 } else { e2 },
            e4a: e2 * e2,
            maxrad: 2.0 * a / f64::EPSILON,
        }
    }

    /// Equatorial radius in meters.
    #[inline]
    pub const fn equatorial_radius(&self) -> f64 {
        self.a
    }

    /// Flattening of the ellipsoid.
    #[inline]
    pub const fn flattening(&self) -> f64 {
        self.f
    }

    /// Convert from geodetic coordinates to geocentric coordinates.
    ///
    /// # Arguments
    ///
    ///  - `lat` - latitude (degrees) [-90., 90.]
    ///  - `lon` - longitude (degrees)
    ///  - `height` - height above the ellipsoid (meters)
    ///
    /// # Returns
    ///
    ///  - `x` - geocentric x (meters).
    ///  - `y` - geocentric y (meters).
    ///  - `z` - geocentric z (meters).
    #[inline]
    pub fn forward(&self, lat: f64, lon: f64, height: f64) -> (f64, f64, f64) {
        let (sin_phi, cos_phi) = geomath::sincosd(geomath::lat_fix(lat));
        let (sin_lam, cos_lam) = geomath::sincosd(lon);
        let n = self.a / (1. - self.e2 * sin_phi * sin_phi).sqrt();
        let r = (n + height) * cos_phi;
        let x = r * cos_lam;
        let y = r * sin_lam;
        let z = (self.e2m * n + height) * sin_phi;
        (x, y, z)
    }

    /// Convert from geodetic coordinates to geocentric coordinates and return a
    /// rotation matrix.
    ///
    /// # Arguments
    ///
    ///  - `lat` - latitude (degrees) [-90., 90.]
    ///  - `lon` - longitude (degrees)
    ///  - `height` - height above the ellipsoid (meters)
    ///
    /// # Returns
    ///
    ///  - `(x, y, z)` - geocentric coordinates (meters).
    ///  - rotation matrix (row-major) mapping a local east, north, up (ENU)
    ///    vector at `(lat, lon, height)` to a geocentric ECEF vector.
    #[inline]
    pub fn forward_with_rotation(
        &self,
        lat: f64,
        lon: f64,
        height: f64,
    ) -> ((f64, f64, f64), [f64; 9]) {
        // Mirrors the `M != nullptr` path of GeographicLib's `IntForward`.
        let (sin_phi, cos_phi) = geomath::sincosd(geomath::lat_fix(lat));
        let (sin_lam, cos_lam) = geomath::sincosd(lon);
        let n = self.a / (1. - self.e2 * sin_phi * sin_phi).sqrt();
        let r = (n + height) * cos_phi;
        let x = r * cos_lam;
        let y = r * sin_lam;
        let z = (self.e2m * n + height) * sin_phi;
        let rotation = rotation_matrix(sin_phi, cos_phi, sin_lam, cos_lam);
        ((x, y, z), rotation)
    }

    /// Convert from geocentric coordinates to geodetic coordinates.
    ///
    /// When multiple geodetic solutions exist, the one minimizing `abs(height)`
    /// is returned.
    ///
    /// The branching mirrors GeographicLib's `Geocentric.cpp`.  The cubic in the
    /// general case uses Vermeille's updated formulas (2011), which are
    /// algebraically equivalent to GeographicLib's `S`/`disc`/`T` staging while
    /// avoiding its `T == 0` guard.
    ///
    /// # Arguments
    ///
    ///  - `x` - geocentric x (meters).
    ///  - `y` - geocentric y (meters).
    ///  - `z` - geocentric z (meters).
    ///
    /// # Returns
    ///
    ///  - `lat` - latitude (degrees).
    ///  - `lon` - longitude (degrees).
    ///  - `height` - height above the ellipsoid (meters).
    #[inline]
    pub fn reverse(&self, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
        let result = self.reverse_internal(x, y, z);
        (result.lat, result.lon, result.height)
    }

    #[inline]
    fn reverse_internal(&self, x: f64, y: f64, z: f64) -> ReverseResult {
        let mut big_r = x.hypot(y);
        let mut slam = if big_r != 0.0 { y / big_r } else { 0.0 };
        let mut clam = if big_r != 0.0 { x / big_r } else { 1.0 };
        let mut h = big_r.hypot(z);
        let (mut sphi, cphi);

        if h > self.maxrad {
            // Match GeographicLib's overflow guard: for extremely distant points,
            // treating the earth as a point gives an adequate height approximation.
            big_r = (x / 2.0).hypot(y / 2.0);
            slam = if big_r != 0.0 { (y / 2.0) / big_r } else { 0.0 };
            clam = if big_r != 0.0 { (x / 2.0) / big_r } else { 1.0 };
            let big_h = (z / 2.0).hypot(big_r);
            sphi = (z / 2.0) / big_h;
            cphi = big_r / big_h;
        } else if self.e4a == 0.0 {
            // Same spherical branch as GeographicLib.  This keeps the origin
            // mapped to the north pole and avoids underflow in the general case.
            let big_h = (if h == 0.0 { 1.0 } else { z }).hypot(big_r);
            sphi = (if h == 0.0 { 1.0 } else { z }) / big_h;
            cphi = big_r / big_h;
            h -= self.a;
        } else {
            let mut p = geomath::sq(big_r / self.a);
            let mut q = self.e2m * geomath::sq(z / self.a);
            let r = (p + q - self.e4a) / 6.0;
            // Match GeographicLib: handle prolate spheroids by swapping the
            // radial and vertical terms before solving the general case.
            if self.f < 0.0 {
                core::mem::swap(&mut p, &mut q);
            }

            if !(self.e4a * q == 0.0 && r <= 0.0) {
                let r3 = r * r * r;
                let e4pq = self.e4a * p * q;
                let evol = 8.0 * r3 + e4pq;

                let u = if evol > 0.0 {
                    // Vermeille (2011) closed form: algebraically equivalent
                    // to GeographicLib's `r + T + r^2/T`, but avoids the sign
                    // pick on T^3 and the `T == 0` special case (the two sqrt
                    // arguments are non-negative and `l > 0` is guaranteed).
                    let l = (evol.sqrt() + e4pq.sqrt()).cbrt();
                    (3.0 * r * r) / (2.0 * l * l) + 0.5 * (l + r / l) * (l + r / l)
                } else {
                    // Vermeille (2011) trigonometric form for the casus
                    // irreducibilis; algebraically equivalent to GeographicLib's
                    // `r + 2*r*cos(ang/3)` branch.
                    let t = 2.0 / 3.0 * e4pq.sqrt().atan2((-evol).sqrt() + (-8.0 * r3).sqrt());
                    -4.0 * r * t.sin() * (FRAC_PI_6 + t).cos()
                };

                let v = (u * u + self.e4a * q).sqrt();
                let uv = if u < 0.0 {
                    self.e4a * q / (v - u)
                } else {
                    u + v
                };
                let w = (self.e2a * (uv - q) / (2.0 * v)).max(0.0);
                let k = uv / ((uv + w * w).sqrt() + w);
                let (k1, k2) = if self.f >= 0.0 {
                    (k, k + self.e2)
                } else {
                    (k - self.e2, k)
                };
                let d = k1 * big_r / k2;
                let big_h = (z / k1).hypot(big_r / k2);
                sphi = (z / k1) / big_h;
                cphi = (big_r / k2) / big_h;
                h = (1.0 - self.e2m / k1) * d.hypot(z);
            } else {
                // Limit branch matching GeographicLib.  The general formula
                // would produce 0/0 for the oblate equatorial plane or the
                // prolate rotation axis.
                let zz = ((if self.f >= 0.0 { self.e4a - p } else { p }) / self.e2m).sqrt();
                let xx = (if self.f < 0.0 { self.e4a - p } else { p }).sqrt();
                let big_h = zz.hypot(xx);
                sphi = zz / big_h;
                cphi = xx / big_h;
                if z < 0.0 {
                    sphi = -sphi;
                }
                h = -self.a * (if self.f >= 0.0 { self.e2m } else { 1.0 }) * big_h / self.e2a;
            }
        }

        ReverseResult {
            lat: geomath::atan2d(sphi, cphi),
            lon: geomath::atan2d(slam, clam),
            height: h,
            sphi,
            cphi,
            slam,
            clam,
        }
    }

    /// Convert from geocentric coordinates to geodetic coordinates and return a
    /// rotation matrix.
    ///
    /// # Arguments
    ///
    ///  - `x` - geocentric x (meters).
    ///  - `y` - geocentric y (meters).
    ///  - `z` - geocentric z (meters).
    ///
    /// # Returns
    ///
    ///  - `(lat, lon, height)` - geodetic coordinates (degrees, degrees, meters).
    ///  - rotation matrix (row-major) mapping a local east, north, up (ENU)
    ///    vector at the returned geodetic position to a geocentric ECEF vector.
    #[inline]
    pub fn reverse_with_rotation(&self, x: f64, y: f64, z: f64) -> ((f64, f64, f64), [f64; 9]) {
        let result = self.reverse_internal(x, y, z);
        let rotation = rotation_matrix(result.sphi, result.cphi, result.slam, result.clam);
        ((result.lat, result.lon, result.height), rotation)
    }
}

#[inline]
fn rotation_matrix(sin_phi: f64, cos_phi: f64, sin_lam: f64, cos_lam: f64) -> [f64; 9] {
    [
        -sin_lam,
        -cos_lam * sin_phi,
        cos_lam * cos_phi,
        cos_lam,
        -sin_lam * sin_phi,
        sin_lam * cos_phi,
        0.,
        cos_phi,
        sin_phi,
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::{assert_abs_diff_eq, assert_relative_eq};
    use std::io::BufRead;

    const GEOCENTRIC_WGS84_TEST_PATH: &str = "test_fixtures/geocentric_wgs84.dat";

    #[test]
    fn sphere() {
        let earth = Geocentric::new(10., 0.);
        assert_eq!(earth.equatorial_radius(), 10.);
        assert_eq!(earth.flattening(), 0.);

        // forward
        {
            let (x, y, z) = earth.forward(0., 0., 5.);
            assert_eq!(x, 15.);
            assert_eq!(y, 0.);
            assert_eq!(z, 0.);
        }

        // reverse
        {
            let (lat, lon, height) = earth.reverse(0., 0., 0.);
            assert_eq!(lat, 90.);
            assert_eq!(lon, 0.);
            assert_eq!(height, -10.);
        }
        {
            let (lat, lon, height) = earth.reverse(10., 0., 0.);
            assert_eq!(lat, 0.);
            assert_eq!(lon, 0.);
            assert_eq!(height, 0.);
        }
    }

    #[test]
    fn reverse_near_origin_with_perturbed_z() {
        let earth = Geocentric::wgs84();

        for (lat, lon, height, z) in [
            (88.10828645_f64, 120.0_f64, -6356728.972246517_f64, 0.0_f64),
            (-88.10828645, 120.0, -6356728.972246517, -f64::MIN_POSITIVE),
        ] {
            let (x, y, _) = earth.forward(lat, lon, height);
            let (lat2, lon2, height2) = earth.reverse(x, y, z);
            assert!((lat - lat2).abs() < 1e-10);
            assert!((lon - lon2).abs() < 1e-10);
            assert!((height - height2).abs() < 30.);
        }
    }

    #[test]
    fn forward_out_of_range_latitude_returns_nan() {
        let earth = Geocentric::wgs84();

        let (x, y, z) = earth.forward(91., 0., 0.);
        assert!(x.is_nan());
        assert!(y.is_nan());
        assert!(z.is_nan());

        let ((x, y, z), rotation) = earth.forward_with_rotation(-91., 0., 0.);
        assert!(x.is_nan());
        assert!(y.is_nan());
        assert!(z.is_nan());
        assert!(rotation.iter().any(|value| value.is_nan()));
    }

    // Mirrors GeographicLib's CartConvert0 / CartConvert1 regression tests
    #[test]
    fn reverse_with_extreme_ellipsoid() {
        // CartConvert 0: heavily oblate (f = 1/100).
        {
            let earth = Geocentric::new(6.4e6, 1. / 100.);
            let (lat, lon, height) = earth.reverse(10e3, 0., 1e3);
            assert!((lat - 85.57).abs() < 0.01);
            assert!((lon - 0.).abs() < 1e-12);
            assert!((height - -6334614.).abs() < 1.);
        }

        // CartConvert 1: prolate (f = -1/100).
        {
            let earth = Geocentric::new(6.4e6, -1. / 100.);
            let (lat, lon, height) = earth.reverse(1e3, 0., 10e3);
            assert!((lat - 4.42).abs() < 0.01);
            assert!((lon - 0.).abs() < 1e-12);
            assert!((height - -6398614.).abs() < 1.);
        }
    }

    #[test]
    fn reverse_far_point_uses_geographiclib_overflow_guard() {
        let earth = Geocentric::wgs84();
        let huge = 1.0e23_f64;

        let (lat, lon, height) = earth.reverse(huge, huge, huge);
        let expected_height = huge.hypot(huge).hypot(huge);

        assert!((lat - 35.264389682754654).abs() < 1e-12);
        assert!((lon - 45.).abs() < 1e-12);
        assert!((height - expected_height).abs() / expected_height < 1e-15);

        let (lat, lon, height) = earth.reverse(f64::MAX, f64::MAX, f64::MAX);

        assert!((lat - 35.264389682754654).abs() < 1e-12);
        assert!((lon - 45.).abs() < 1e-12);
        assert!(height.is_infinite());

        let (lat, lon, height) = earth.reverse(0., 0., huge);

        assert_eq!(lat, 90.);
        assert_eq!(lon, 0.);
        assert_eq!(height, huge);
    }

    #[test]
    fn prolate_roundtrip() {
        let earth = Geocentric::new(6.4e6, -1. / 100.);
        let (lat, lon, height) = (35., 123., 1000.);

        let (x, y, z) = earth.forward(lat, lon, height);
        let (lat2, lon2, height2) = earth.reverse(x, y, z);

        assert!((lat - lat2).abs() < 1e-10);
        assert!((lon - lon2).abs() < 1e-10);
        assert!((height - height2).abs() < 1e-7);

        let (lat, lon, height) = earth.reverse(0., 0., 0.);
        assert_eq!(lat, 0.);
        assert_eq!(lon, 0.);
        assert!((height - -6.4e6).abs() < 1.0);
    }

    #[test]
    fn forward_with_rotation() {
        let earth = Geocentric::wgs84();
        let (_xyz, m) = earth.forward_with_rotation(0., 0., 0.);
        assert_eq!(m, [-0., -0., 1., 1., -0., 0., 0., 1., 0.]);
    }

    #[test]
    fn reverse_with_rotation_matches_forward_rotation() {
        let earth = Geocentric::wgs84();
        let ((x, y, z), forward_rotation) = earth.forward_with_rotation(37., 140., 50.);
        let (_lla, reverse_rotation) = earth.reverse_with_rotation(x, y, z);

        for (forward, reverse) in forward_rotation.iter().zip(reverse_rotation.iter()) {
            assert!((forward - reverse).abs() < 1e-12);
        }
    }

    #[test]
    fn geocentric_cpp_compatibility() {
        // Generated from C++ GeographicLib Geocentric with std::setprecision(17).
        let file = std::fs::File::open(GEOCENTRIC_WGS84_TEST_PATH)
            .expect("failed to open geocentric_wgs84.dat");
        let earth = Geocentric::wgs84();
        let reader = std::io::BufReader::new(file);

        for (i, line) in reader.lines().enumerate() {
            let line_num = i + 1;
            let line = line.expect("failed to read GeocentricTest line");
            let values: Vec<f64> = line
                .split_whitespace()
                .map(|value| {
                    value
                        .parse::<f64>()
                        .expect("failed to parse GeocentricTest value")
                })
                .collect();
            assert_eq!(values.len(), 9, "line {line_num}");

            let (x, y, z) = earth.forward(values[0], values[1], values[2]);
            assert_relative_eq!(x, values[3], epsilon = 1e-8, max_relative = 1e-14);
            assert_relative_eq!(y, values[4], epsilon = 1e-8, max_relative = 1e-14);
            assert_relative_eq!(z, values[5], epsilon = 1e-8, max_relative = 1e-14);

            let (lat, lon, height) = earth.reverse(values[3], values[4], values[5]);
            assert_abs_diff_eq!(lat, values[6], epsilon = 1e-12);
            assert_abs_diff_eq!(lon, values[7], epsilon = 1e-12);
            assert_relative_eq!(height, values[8], epsilon = 1e-7, max_relative = 1e-14);
        }
    }
}
