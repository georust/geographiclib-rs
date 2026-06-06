//! Transverse Mercator projection.

use crate::geodesic::{WGS84_A, WGS84_F};
use crate::geomath;
use std::f64::consts::PI;

const TRANSVERSE_MERCATOR_ORDER: usize = 6;
const UTM_K0: f64 = 0.9996;

type Coeffs = [f64; TRANSVERSE_MERCATOR_ORDER];

const B1_COEFF: [f64; TRANSVERSE_MERCATOR_ORDER / 2 + 2] = [1.0, 4.0, 64.0, 256.0, 256.0];

const ALP_COEFF: [f64; TRANSVERSE_MERCATOR_ORDER * (TRANSVERSE_MERCATOR_ORDER + 3) / 2] = [
    31564.0,
    -66675.0,
    34440.0,
    47250.0,
    -100800.0,
    75600.0,
    151200.0,
    -1983433.0,
    863232.0,
    748608.0,
    -1161216.0,
    524160.0,
    1935360.0,
    670412.0,
    406647.0,
    -533952.0,
    184464.0,
    725760.0,
    6601661.0,
    -7732800.0,
    2230245.0,
    7257600.0,
    -13675556.0,
    3438171.0,
    7983360.0,
    212378941.0,
    319334400.0,
];

const BET_COEFF: [f64; TRANSVERSE_MERCATOR_ORDER * (TRANSVERSE_MERCATOR_ORDER + 3) / 2] = [
    384796.0,
    -382725.0,
    -6720.0,
    932400.0,
    -1612800.0,
    1209600.0,
    2419200.0,
    -1118711.0,
    1695744.0,
    -1174656.0,
    258048.0,
    80640.0,
    3870720.0,
    22276.0,
    -16929.0,
    -15984.0,
    12852.0,
    362880.0,
    -830251.0,
    -158400.0,
    197865.0,
    7257600.0,
    -435388.0,
    453717.0,
    15966720.0,
    20648693.0,
    638668800.0,
];

/// Transverse Mercator projection.
///
/// This uses Krüger's method which evaluates the
/// projection and its inverse in terms of a series.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct TransverseMercator {
    a: f64,
    f: f64,
    k0: f64,
    e2: f64,
    e2m: f64,
    es: f64,
    c: f64,
    a1: f64,
    b1: f64,
    bet: Coeffs,
    alp: Coeffs,
}

impl TransverseMercator {
    /// Creates a new Transverse Mercator projection.
    ///
    /// # Arguments
    ///
    /// - `a`: equatorial radius in meters.
    /// - `f`: flattening of the ellipsoid.
    /// - `k0`: central scale factor.
    pub fn new(a: f64, f: f64, k0: f64) -> Self {
        let e2 = f * (2.0 - f);
        let e2m = 1.0 - e2;
        let es = if f < 0.0 { -1.0 } else { 1.0 } * e2.abs().sqrt();
        let c = e2m.sqrt() * geomath::eatanhe(1.0, es).exp();
        let (b1, bet, alp) = setup_coefficients(f);
        let a1 = a * b1;
        TransverseMercator {
            a,
            f,
            k0,
            e2,
            e2m,
            es,
            c,
            a1,
            b1,
            bet,
            alp,
        }
    }

    /// Creates the UTM Transverse Mercator projection for WGS84.
    ///
    /// This does not add UTM false easting or false northing.
    pub fn utm() -> Self {
        TransverseMercator::new(WGS84_A, WGS84_F, UTM_K0)
    }

    /// Returns the equatorial radius in meters.
    pub fn equatorial_radius(&self) -> f64 {
        self.a
    }

    /// Returns the ellipsoid flattening.
    pub fn flattening(&self) -> f64 {
        self.f
    }

    /// Returns the central scale factor.
    pub fn central_scale(&self) -> f64 {
        self.k0
    }

    /// Forward projection. Returns `(x, y, gamma, k)`.
    pub fn forward(&self, lon0: f64, lat: f64, lon: f64) -> (f64, f64, f64, f64) {
        let mut lat = geomath::lat_fix(lat);
        let (mut lon, _lon_error) = geomath::ang_diff(lon0, lon);

        let mut latsign: f64 = 1.0_f64.copysign(lat);
        let lonsign: f64 = 1.0_f64.copysign(lon);
        lat *= latsign;
        lon *= lonsign;

        let backside = lon > 90.0;
        if backside {
            if lat == 0.0 {
                latsign = -1.0;
            }
            lon = 180.0 - lon;
        }

        let at_pole = lat == 90.0;
        let (eta, xi, mut gamma, mut k) = self.forward_inner(lat, lon, at_pole);
        let scale = self.a1 * self.k0;
        let x = scale * eta * lonsign;
        let y = scale * (if backside { PI - xi } else { xi }) * latsign;

        if backside {
            gamma = 180.0 - gamma;
        }
        gamma *= latsign * lonsign;
        gamma = geomath::ang_normalize(gamma);
        k *= self.k0;

        (x, y, gamma, k)
    }

    /// Reverse projection. Returns `(lat, lon, gamma, k)`.
    pub fn reverse(&self, lon0: f64, x: f64, y: f64) -> (f64, f64, f64, f64) {
        let denom = self.a1 * self.k0;
        let mut xi = y / denom;
        let mut eta = x / denom;

        let xisign: f64 = 1.0_f64.copysign(xi);
        let etasign: f64 = 1.0_f64.copysign(eta);
        xi *= xisign;
        eta *= etasign;

        let backside = xi > PI / 2.0;
        if backside {
            xi = PI - xi;
        }

        let (mut lat, mut lon, mut gamma, mut k) = self.reverse_inner(xi, eta);
        lat *= xisign;
        if backside {
            lon = PI - lon;
        }
        lon *= etasign;
        if backside {
            gamma = 180.0 - gamma;
        }
        gamma *= xisign * etasign;
        gamma = geomath::ang_normalize(gamma);
        k *= self.k0;

        (
            lat.to_degrees(),
            geomath::ang_normalize(lon.to_degrees() + lon0),
            gamma,
            k,
        )
    }

    fn forward_inner(&self, lat: f64, lon: f64, at_pole: bool) -> (f64, f64, f64, f64) {
        let (xip, etap, mut gamma, mut k) = if at_pole {
            (PI / 2.0, 0.0, lon, self.c)
        } else {
            let (sphi, cphi) = geomath::sincosd(lat);
            let (slam, clam) = geomath::sincosd(lon);
            let tau = sphi / cphi;
            let taup = geomath::taupf(tau, self.es);
            let xip = f64::atan2(taup, clam);
            let etap = (slam / f64::hypot(taup, clam)).asinh();
            // Krueger p 22 (44)
            let gamma = geomath::atan2d(slam * taup, clam * f64::hypot(1.0, taup));
            // k0 = sqrt(1 - e2 * sin(phi)^2) * (cos(phi') / cos(phi)) * cosh(etap)
            let k = (self.e2m + self.e2 * cphi * cphi).sqrt() * f64::hypot(1.0, tau)
                / f64::hypot(taup, clam);
            (xip, etap, gamma, k)
        };

        let (s0, c0) = (2.0 * xip).sin_cos();
        let sh0 = (2.0 * etap).sinh();
        let ch0 = (2.0 * etap).cosh();

        let ((dxi, deta), (dz_re, dz_im)) = clenshaw(&self.alp, s0, c0, sh0, ch0);
        let xi = xip + dxi;
        let eta = etap + deta;

        gamma -= geomath::atan2d(dz_im, dz_re);
        k *= self.b1 * f64::hypot(dz_re, dz_im);

        (eta, xi, gamma, k)
    }

    fn reverse_inner(&self, xi: f64, eta: f64) -> (f64, f64, f64, f64) {
        let (s0, c0) = (2.0 * xi).sin_cos();
        let sh0 = (2.0 * eta).sinh();
        let ch0 = (2.0 * eta).cosh();

        let ((dxi, deta), (dz_re, dz_im)) = clenshaw(&self.bet, s0, c0, sh0, ch0);
        let xip = xi + dxi;
        let etap = eta + deta;

        let mut gamma = geomath::atan2d(dz_im, dz_re);
        let mut k = self.b1 / f64::hypot(dz_re, dz_im);

        let s = etap.sinh();
        // cos(pi/2) might be negative
        let c = f64::max(0.0, xip.cos());
        let r = f64::hypot(s, c);

        let (lat, lon) = if r != 0.0 {
            let lon = f64::atan2(s, c);
            let sxip = xip.sin();
            let tau = geomath::tauf(sxip / r, self.es);
            // Krueger p 19 (31)
            gamma += geomath::atan2d(sxip * etap.tanh(), c);
            // Note cos(phi') * cosh(eta') = r
            k *= (self.e2m + self.e2 / (1.0 + tau * tau)).sqrt() * f64::hypot(1.0, tau) * r;
            (tau.atan(), lon)
        } else {
            k *= self.c;
            (PI / 2.0, 0.0)
        };
        (lat, lon, gamma, k)
    }
}

fn setup_coefficients(f: f64) -> (f64, Coeffs, Coeffs) {
    let n = f / (2.0 - f);

    let mut bet = [0.0; TRANSVERSE_MERCATOR_ORDER];
    let mut alp = [0.0; TRANSVERSE_MERCATOR_ORDER];

    let m = TRANSVERSE_MERCATOR_ORDER / 2;
    let b1 = geomath::polyval(m, &B1_COEFF, geomath::sq(n)) / (B1_COEFF[m + 1] * (1.0 + n));

    let mut offset = 0;
    let mut d = n;
    #[allow(clippy::needless_range_loop)]
    for l in 0..TRANSVERSE_MERCATOR_ORDER {
        let m = TRANSVERSE_MERCATOR_ORDER - l - 1;
        alp[l] = d * geomath::polyval(m, &ALP_COEFF[offset..], n) / ALP_COEFF[offset + m + 1];
        bet[l] = -d * geomath::polyval(m, &BET_COEFF[offset..], n) / BET_COEFF[offset + m + 1];
        offset += m + 2;
        d *= n;
    }

    (b1, bet, alp)
}

/// Clenshaw summation of the Krüger series and its derivative.
fn clenshaw(a: &Coeffs, s0: f64, c0: f64, sh0: f64, ch0: f64) -> ((f64, f64), (f64, f64)) {
    let n = a.len();

    // (r2, i2) is Karney's `a = 2 * cos(2*zeta)` recurrence multiplier.
    let r2 = 2.0 * c0 * ch0;
    let i2 = -2.0 * s0 * sh0;

    let (mut y_re, mut y_im) = (a[n - 1], 0.0);
    let (mut y1_re, mut y1_im) = (0.0, 0.0);
    let (mut z_re, mut z_im) = ((2 * n) as f64 * a[n - 1], 0.0);
    let (mut z1_re, mut z1_im) = (0.0, 0.0);

    for (idx, v) in a[..n - 1].iter().enumerate().rev() {
        let zv = 2.0 * (idx + 1) as f64 * v;

        let (y2_re, y2_im) = (y1_re, y1_im);
        (y1_re, y1_im) = (y_re, y_im);
        y_re = -y2_re + r2 * y1_re - i2 * y1_im + v;
        y_im = -y2_im + i2 * y1_re + r2 * y1_im;

        let (z2_re, z2_im) = (z1_re, z1_im);
        (z1_re, z1_im) = (z_re, z_im);
        z_re = -z2_re + r2 * z1_re - i2 * z1_im + zv;
        z_im = -z2_im + i2 * z1_re + r2 * z1_im;
    }

    // sin(2*zeta)
    let sin_re = s0 * ch0;
    let sin_im = c0 * sh0;
    let sum_re = sin_re * y_re - sin_im * y_im;
    let sum_im = sin_re * y_im + sin_im * y_re;

    // cos(2*zeta)
    let cos_re = c0 * ch0;
    let cos_im = -s0 * sh0;
    let deriv_re = 1.0 - z1_re + cos_re * z_re - cos_im * z_im;
    let deriv_im = -z1_im + cos_re * z_im + cos_im * z_re;

    ((sum_re, sum_im), (deriv_re, deriv_im))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    // Reference values generated by GeographicLib (C++) printed with %.17g.

    #[test]
    fn round_trip() {
        let tm = TransverseMercator::utm();
        assert_eq!(
            (tm.equatorial_radius(), tm.flattening(), tm.central_scale()),
            (WGS84_A, WGS84_F, UTM_K0),
        );
        let (x, y, gamma, k) = tm.forward(139.0, 35.0, 139.25);

        assert_relative_eq!(x, 22812.938572148290, max_relative = 1e-15);
        assert_relative_eq!(y, 3873071.6116043567, max_relative = 1e-15);
        assert_relative_eq!(gamma, 0.14339472802333403, max_relative = 1e-15);
        assert_relative_eq!(k, 0.99960641388268068, max_relative = 1e-15);

        let (lat, lon, gamma_r, k_r) = tm.reverse(139.0, x, y);
        assert_relative_eq!(lat, 35.0, epsilon = 1e-13);
        assert_relative_eq!(lon, 139.25, max_relative = 1e-15);
        assert_relative_eq!(gamma_r, gamma, max_relative = 1e-13);
        assert_relative_eq!(k_r, k, max_relative = 1e-15);

        // Southern hemisphere: exercises latsign=-1 (forward) and xisign=-1 (reverse).
        let (xs, ys, gammas, ks) = tm.forward(139.0, -35.0, 139.25);
        assert_relative_eq!(xs, x, max_relative = 1e-15);
        assert_relative_eq!(ys, -y, max_relative = 1e-15);
        assert_relative_eq!(gammas, -gamma, max_relative = 1e-15);
        assert_relative_eq!(ks, k, max_relative = 1e-15);

        let (lats, lons, _, _) = tm.reverse(139.0, xs, ys);
        assert_relative_eq!(lats, -35.0, epsilon = 1e-13);
        assert_relative_eq!(lons, 139.25, max_relative = 1e-15);
    }

    // GeographicLib supports the backside of the projection beyond 90 deg from lon0.
    #[test]
    fn round_trip_on_backside() {
        let tm = TransverseMercator::utm();
        let (x, y, gamma, k) = tm.forward(0.0, 10.0, 179.0);

        assert_relative_eq!(x, 109600.77251445431, max_relative = 1e-15);
        assert_relative_eq!(y, 18890351.296849594, max_relative = 1e-15);
        assert_relative_eq!(gamma, 179.8263343830516, max_relative = 1e-15);
        assert_relative_eq!(k, 0.99974864012035936, max_relative = 1e-15);

        let (lat, lon, _, _) = tm.reverse(0.0, x, y);
        assert_relative_eq!(lat, 10.0, epsilon = 1e-13);
        assert_relative_eq!(lon, 179.0, max_relative = 1e-15);
    }

    #[test]
    fn forward_on_equator_backside() {
        let tm = TransverseMercator::utm();
        let (x, y, gamma, k) = tm.forward(0.0, 0.0, 179.0);

        assert_relative_eq!(x, 111280.65089140118, max_relative = 1e-15);
        assert_relative_eq!(y, -19995929.886041995, max_relative = 1e-15);
        assert_relative_eq!(gamma, -180.0, max_relative = 1e-15);
        assert_relative_eq!(k, 0.99975329355317955, max_relative = 1e-15);
    }

    #[test]
    fn reverse_spherical_pole() {
        let tm = TransverseMercator::new(WGS84_A, 0.0, 1.0);
        let (lat, lon, gamma, k) = tm.reverse(0.0, 0.0, WGS84_A * PI / 2.0);

        assert_relative_eq!(lat, 89.999999999999986, max_relative = 1e-15);
        assert_relative_eq!(lon, 180.0, max_relative = 1e-15);
        assert_relative_eq!(gamma, 180.0, max_relative = 1e-15);
        assert_relative_eq!(k, 1.0, max_relative = 1e-15);
    }

    #[test]
    fn reverse_inner_handles_rounded_pole() {
        let tm = TransverseMercator::new(WGS84_A, 0.0, 1.0);
        let xi = f64::from_bits((PI / 2.0).to_bits() + 1);
        let (lat_rad, lon_rad, _, k) = tm.reverse_inner(xi, 0.0);

        assert_relative_eq!(lat_rad, PI / 2.0, max_relative = 1e-15);
        assert_relative_eq!(lon_rad, 0.0, max_relative = 1e-15);
        assert_relative_eq!(k, 1.0, max_relative = 1e-15);
    }

    #[test]
    fn prolate_ellipsoid_round_trip() {
        let tm = TransverseMercator::new(6.4e6, -1.0 / 298.257223563, 1.0);
        let (x, y, _, _) = tm.forward(0.0, 45.0, 1.0);
        let (lat, lon, _, _) = tm.reverse(0.0, x, y);
        assert_relative_eq!(lat, 45.0, epsilon = 1e-11);
        assert_relative_eq!(lon, 1.0, max_relative = 1e-15);
    }

    // From GeographicLib TransverseMercatorProj1: meridian convergence regression at the pole (2013-06-26).
    #[test]
    fn geographiclib_proj1_forward_at_pole() {
        let tm = TransverseMercator::new(WGS84_A, WGS84_F, 1.0);
        let (x, y, gamma, k) = tm.forward(0.0, 90.0, 75.0);
        assert_eq!(x, 0.0); // exact: at_pole branch produces eta=0
        assert_relative_eq!(y, 10001965.729312722, max_relative = 1e-15);
        assert_relative_eq!(gamma, 75.0, max_relative = 1e-15);
        assert_relative_eq!(k, 0.99999999999999978, max_relative = 1e-15);
    }

    // From GeographicLib TransverseMercatorProj3: scale regression at the pole (2013-06-30).
    #[test]
    fn geographiclib_proj3_reverse_at_quarter_meridian() {
        let tm = TransverseMercator::new(WGS84_A, WGS84_F, 1.0);
        let (lat, lon, gamma, k) = tm.reverse(0.0, 0.0, 10001965.7293127228);
        assert_relative_eq!(lat, 89.999999999999972, max_relative = 1e-15);
        assert_relative_eq!(lon, 180.0, max_relative = 1e-15);
        assert_relative_eq!(gamma, 180.0, max_relative = 1e-15);
        assert_relative_eq!(k, 0.99999999999999967, max_relative = 1e-15);
    }

    // From GeographicLib TransverseMercatorProj5: complex Clenshaw on a heavily oblate ellipsoid (2017-04-15).
    #[test]
    fn geographiclib_proj5_forward_high_flattening() {
        let tm = TransverseMercator::new(6.4e6, 1.0 / 150.0, 0.9996);
        let (x, y, gamma, k) = tm.forward(0.0, 20.0, 30.0);
        assert_relative_eq!(x, 3266035.4538597651, max_relative = 1e-15);
        assert_relative_eq!(y, 2518371.5526758581, max_relative = 1e-15);
        assert_relative_eq!(gamma, 11.20735650214103, max_relative = 1e-15);
        assert_relative_eq!(k, 1.1341389607409129, max_relative = 1e-15);
    }

    // From GeographicLib TransverseMercatorProj7: reverse counterpart on the same heavily oblate ellipsoid.
    #[test]
    fn geographiclib_proj7_reverse_high_flattening() {
        let tm = TransverseMercator::new(6.4e6, 1.0 / 150.0, 0.9996);
        let (lat, lon, gamma, k) = tm.reverse(0.0, 3.3e6, 2.5e6);
        assert_relative_eq!(lat, 19.80370996792691, max_relative = 1e-15);
        assert_relative_eq!(lon, 30.249197022823417, max_relative = 1e-15);
        assert_relative_eq!(gamma, 11.214378172893015, max_relative = 1e-15);
        assert_relative_eq!(k, 1.1370257757585773, max_relative = 1e-15);
    }

    // Mirrors GeographicLib's develop/NaNTester.cpp.
    #[test]
    fn nan_propagation() {
        let tm = TransverseMercator::utm();
        let nan = f64::NAN;

        for (x, y, gamma, k) in [
            tm.forward(nan, 0.0, 0.0),
            tm.forward(0.0, nan, 0.0),
            tm.forward(0.0, 0.0, nan),
        ] {
            assert!(x.is_nan() && y.is_nan() && gamma.is_nan() && k.is_nan());
        }

        for (lat, lon, gamma, k) in [tm.reverse(0.0, nan, 0.0), tm.reverse(0.0, 0.0, nan)] {
            assert!(lat.is_nan() && lon.is_nan() && gamma.is_nan() && k.is_nan());
        }

        // NaN lon0 only contaminates lon (matches GeographicLib).
        let (lat, lon, gamma, k) = tm.reverse(nan, 0.0, 0.0);
        assert!(!lat.is_nan() && lon.is_nan() && !gamma.is_nan() && !k.is_nan());
    }

    // Karney's TMcoords.dat fixture. See test_fixtures/README.md.
    // Rows requiring `TransverseMercatorExact` are skipped.
    #[test]
    fn tmcoords_compatibility() {
        use std::io::BufRead;

        const BUILTIN_TEST_PATH: &str = "test_fixtures/TMcoords-excerpt.dat";
        const FULL_TEST_PATH: &str = "test_fixtures/test_data_unzipped/TMcoords.dat";

        let full = cfg!(feature = "test_full");
        let path = if full {
            FULL_TEST_PATH
        } else {
            BUILTIN_TEST_PATH
        };

        // Karney's TMcoords.dat uses WGS84, central meridian 0, k0 = 0.9996.
        let tm = TransverseMercator::utm();
        let file = std::fs::File::open(path).unwrap_or_else(|_| {
            panic!(
                "failed to open {}; run script/download-test-data.sh for the full fixture",
                path
            )
        });
        let reader = std::io::BufReader::new(file);

        for (i, line) in reader.lines().enumerate() {
            let line_num = i + 1;
            let line = line.expect("failed to read TMcoords line");
            let v: Vec<f64> = line
                .split_whitespace()
                .map(|s| s.parse::<f64>().expect("failed to parse TMcoords value"))
                .collect();
            assert_eq!(v.len(), 6, "line {line_num}");
            let (lat, lon, x_ref, y_ref, gamma_ref, k_ref) = (v[0], v[1], v[2], v[3], v[4], v[5]);

            // For the full fixture, skip rows from segments we cannot resolve.
            #[cfg(feature = "test_full")]
            if !krueger_compatible(line_num, lat, lon) {
                continue;
            }
            #[cfg(not(feature = "test_full"))]
            let _ = line_num;

            let (x, y, gamma, k) = tm.forward(0.0, lat, lon);
            assert_relative_eq!(x, x_ref, epsilon = 1e-9, max_relative = 1e-13);
            assert_relative_eq!(y, y_ref, epsilon = 1e-9, max_relative = 1e-13);
            assert_relative_eq!(gamma, gamma_ref, epsilon = 1e-14, max_relative = 1e-14);
            assert_relative_eq!(k, k_ref, max_relative = 1e-15);
        }
    }

    /// Returns false for `TMcoords.dat` rows that require `TransverseMercatorExact`
    /// (not implemented here) and cannot be resolved by the Krueger 6th-order series.
    #[cfg(feature = "test_full")]
    fn krueger_compatible(line_num: usize, lat: f64, lon: f64) -> bool {
        // Segments 7..13 (rows 255001..=287000) cover the Krueger branch cut
        // and the lat < 0 "extended" domain.
        if line_num > 255000 {
            return false;
        }
        // Mirror ETA_MAX in script/excerpt-tmcoords.py.
        let (sphi, cphi) = lat.to_radians().sin_cos();
        let (slam, clam) = lon.to_radians().sin_cos();
        let h = f64::hypot(sphi / cphi, clam);
        (slam.abs() / h).asinh() < 0.5
    }
}
