#![allow(non_snake_case)]

use std::cmp::Ordering;

use crate::geodesiccapability::{
    AZIMUTH, DISTANCE_IN, GEODESICSCALE, LATITUDE, LONGITUDE, REDUCEDLENGTH,
};
use crate::geodesicline::GeodesicLine;
use crate::geomath;
use crate::Geodesic;

const NUMIT_: usize = 20;

pub struct Gnomonic {
    eps0_: f64,
    eps_: f64,
    _earth: Geodesic,
    _a: f64,
    _f: f64,
}

impl Gnomonic {
    pub fn wgs84() -> Self {
        let _earth = Geodesic::wgs84();
        let eps0_ = f64::EPSILON;
        Gnomonic {
            eps0_,
            eps_: 0.01 * eps0_.sqrt(),
            _earth,
            _a: _earth.equatorial_radius(),
            _f: _earth.flattening(),
        }
    }

    pub fn new(eps0_: f64, eps_: f64, _earth: Geodesic, _a: f64, _f: f64) -> Self {
        Gnomonic {
            eps0_,
            eps_,
            _earth,
            _a,
            _f,
        }
    }

    pub fn equatorial_radius(&self) -> f64 {
        self._earth.equatorial_radius()
    }

    pub fn flattening(&self) -> f64 {
        self._earth.flattening()
    }
}

impl Gnomonic {
    // returns (x, y, azi, rk)
    pub fn forward(&self, lat0: f64, lon0: f64, lat: f64, lon: f64) -> (f64, f64, f64, f64) {
        let mask = AZIMUTH | REDUCEDLENGTH | GEODESICSCALE;
        let (_, _, azi0, azi, m, M, _, _) =
            self._earth._gen_inverse_azi(lat0, lon0, lat, lon, mask);

        let rk = M;
        let (x, y) = if M <= 0.0 {
            (f64::NAN, f64::NAN)
        } else {
            let rho = m / M;
            let (x, y) = geomath::sincosd(azi0);
            (x * rho, y * rho)
        };

        (x, y, azi, rk)
    }

    // returns (lat, lon, azi, rk)
    pub fn reverse(&self, lat0: f64, lon0: f64, x: f64, y: f64) -> (f64, f64, f64, f64) {
        let azi0 = geomath::atan2d(x, y);
        let mut rho = f64::hypot(x, y);
        let mut s = self._a * (rho / self._a).atan();
        let little = rho <= self._a;

        if !little {
            rho = 1.0 / rho;
        }

        let mask = LATITUDE | LONGITUDE | AZIMUTH | DISTANCE_IN | REDUCEDLENGTH | GEODESICSCALE;
        let line = GeodesicLine::new(&self._earth, lat0, lon0, azi0, Some(mask), None, None);

        let mut trip = 0;
        let mut lat1 = f64::NAN;
        let mut lon1 = f64::NAN;
        let mut azi1 = f64::NAN;
        let mut M = f64::NAN; // rk

        for _ in 0..NUMIT_ {
            let mask = LATITUDE | LONGITUDE | AZIMUTH | REDUCEDLENGTH | GEODESICSCALE;
            let result = line._gen_position(false, s, mask);
            let m = result.5;
            lat1 = result.1;
            lon1 = result.2;
            azi1 = result.3;
            M = result.6;

            if trip > 0 {
                break;
            }

            // If little, solve rho(s) = rho with drho(s)/ds = 1/M^2
            // else solve 1/rho(s) = 1/rho with d(1/rho(s))/ds = -1/m^2
            let ds = if little {
                (m - rho * M) * M
            } else {
                (rho * m - M) * m
            };

            s -= ds;
            // if !(ds.abs() >= self.eps_ * self._a) {
            if let None | Some(Ordering::Less) = ds.abs().partial_cmp(&(self.eps_ * self._a)) {
                trip += 1;
            }
        }

        (lat1, lon1, azi1, M)
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn forward() {
        // Paris
        let lat0 = 48.0 + 50.0 / 60.0;
        let lon0 = 2.0 + 20.0 / 60.0;

        let proj = Gnomonic::wgs84();

        // Calais
        let lat = 50.9;
        let lon = 1.8;

        let (x, y, _, _) = proj.forward(lat0, lon0, lat, lon);

        assert_eq!(x, -37543.66988338346);
        assert_eq!(y, 230103.2138503046);
    }

    #[test]
    fn reverse() {
        // Paris
        let lat0 = 48.0 + 50.0 / 60.0;
        let lon0 = 2.0 + 20.0 / 60.0;

        let proj = Gnomonic::wgs84();

        let x = -38e3;
        let y = 230e3;

        let (lat, lon, _, _) = proj.reverse(lat0, lon0, x, y);

        assert_eq!(lat, 50.899043656956295);
        assert_eq!(lon, 1.7935283011976253);
    }
}
