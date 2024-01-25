#![allow(non_snake_case)]

use crate::geodesic::{self, GEODESIC_ORDER};
use crate::geodesic_capability as caps;
use crate::geomath;
use std::collections::HashMap;

#[derive(Copy, Clone, PartialEq, PartialOrd, Debug)]
pub struct GeodesicLine {
    tiny_: f64, // This should be moved to consts
    _A1m1: f64,
    _A2m1: f64,
    _A3c: f64,
    _A4: f64,
    _B11: f64,
    _B21: f64,
    _B31: f64,
    _B41: f64,
    _C1a: [f64; GEODESIC_ORDER + 1],
    _C1pa: [f64; GEODESIC_ORDER + 1],
    _C2a: [f64; GEODESIC_ORDER + 1],
    _C3a: [f64; GEODESIC_ORDER],
    _C4a: [f64; GEODESIC_ORDER],
    _b: f64,
    _c2: f64,
    _calp0: f64,
    _csig1: f64,
    _comg1: f64,
    _ctau1: f64,
    _dn1: f64,
    _f1: f64,
    _k2: f64,
    _salp0: f64,
    _somg1: f64,
    _ssig1: f64,
    _stau1: f64,
    _a13: f64,
    _a: f64,
    azi1: f64,
    calp1: f64,
    caps: u64,
    f: f64,
    lat1: f64,
    lon1: f64,
    _s13: f64,
    salp1: f64,
}

impl GeodesicLine {
    pub fn new(
        geod: &geodesic::Geodesic,
        lat1: f64,
        lon1: f64,
        azi1: f64,
        caps: Option<u64>,
        salp1: Option<f64>,
        calp1: Option<f64>,
    ) -> Self {
        let caps = match caps {
            None => caps::STANDARD | caps::DISTANCE_IN,
            Some(caps) => caps,
        };
        let salp1 = match salp1 {
            None => std::f64::NAN,
            Some(salp1) => salp1,
        };
        let calp1 = match calp1 {
            None => std::f64::NAN,
            Some(calp1) => calp1,
        };

        // This was taken from geodesic, putting it here for convenience
        let tiny_ = geomath::get_min_val().sqrt();

        let _a = geod.a;
        let f = geod.f;
        let _b = geod._b;
        let _c2 = geod._c2;
        let _f1 = geod._f1;
        let caps = caps | caps::LATITUDE | caps::AZIMUTH | caps::LONG_UNROLL;
        let (azi1, salp1, calp1) = if salp1.is_nan() || calp1.is_nan() {
            let azi1 = geomath::ang_normalize(azi1);
            let (salp1, calp1) = geomath::sincosd(geomath::ang_round(azi1));
            (azi1, salp1, calp1)
        } else {
            (azi1, salp1, calp1)
        };
        let lat1 = geomath::lat_fix(lat1);

        let (mut sbet1, mut cbet1) = geomath::sincosd(geomath::ang_round(lat1));
        sbet1 *= _f1;
        geomath::norm(&mut sbet1, &mut cbet1);
        cbet1 = tiny_.max(cbet1);
        let _dn1 = (1.0 + geod._ep2 * geomath::sq(sbet1)).sqrt();
        let _salp0 = salp1 * cbet1;
        let _calp0 = calp1.hypot(salp1 * sbet1);
        let mut _ssig1 = sbet1;
        let _somg1 = _salp0 * sbet1;
        let mut _csig1 = if sbet1 != 0.0 || calp1 != 0.0 {
            cbet1 * calp1
        } else {
            1.0
        };
        let _comg1 = _csig1;
        geomath::norm(&mut _ssig1, &mut _csig1);
        let _k2 = geomath::sq(_calp0) * geod._ep2;
        let eps = _k2 / (2.0 * (1.0 + (1.0 + _k2).sqrt()) + _k2);

        let mut _A1m1 = 0.0;
        let mut _C1a: [f64; GEODESIC_ORDER + 1] = [0.0; GEODESIC_ORDER + 1];
        let mut _B11 = 0.0;
        let mut _stau1 = 0.0;
        let mut _ctau1 = 0.0;
        if caps & caps::CAP_C1 != 0 {
            _A1m1 = geomath::_A1m1f(eps, geod.GEODESIC_ORDER);
            geomath::_C1f(eps, &mut _C1a, geod.GEODESIC_ORDER);
            _B11 = geomath::sin_cos_series(true, _ssig1, _csig1, &_C1a);
            let s = _B11.sin();
            let c = _B11.cos();
            _stau1 = _ssig1 * c + _csig1 * s;
            _ctau1 = _csig1 * c - _ssig1 * s;
        }

        let mut _C1pa: [f64; GEODESIC_ORDER + 1] = [0.0; GEODESIC_ORDER + 1];
        if caps & caps::CAP_C1p != 0 {
            geomath::_C1pf(eps, &mut _C1pa, geod.GEODESIC_ORDER);
        }

        let mut _A2m1 = 0.0;
        let mut _C2a: [f64; GEODESIC_ORDER + 1] = [0.0; GEODESIC_ORDER + 1];
        let mut _B21 = 0.0;
        if caps & caps::CAP_C2 != 0 {
            _A2m1 = geomath::_A2m1f(eps, geod.GEODESIC_ORDER);
            geomath::_C2f(eps, &mut _C2a, geod.GEODESIC_ORDER);
            _B21 = geomath::sin_cos_series(true, _ssig1, _csig1, &_C2a);
        }

        let mut _C3a: [f64; GEODESIC_ORDER] = [0.0; GEODESIC_ORDER];
        let mut _A3c = 0.0;
        let mut _B31 = 0.0;
        if caps & caps::CAP_C3 != 0 {
            geod._C3f(eps, &mut _C3a);
            _A3c = -f * _salp0 * geod._A3f(eps);
            _B31 = geomath::sin_cos_series(true, _ssig1, _csig1, &_C3a);
        }

        let mut _C4a: [f64; GEODESIC_ORDER] = [0.0; GEODESIC_ORDER];
        let mut _A4 = 0.0;
        let mut _B41 = 0.0;
        if caps & caps::CAP_C4 != 0 {
            geod._C4f(eps, &mut _C4a);
            _A4 = geomath::sq(_a) * _calp0 * _salp0 * geod._e2;
            _B41 = geomath::sin_cos_series(false, _ssig1, _csig1, &_C4a);
        }

        let _s13 = std::f64::NAN;
        let _a13 = std::f64::NAN;

        GeodesicLine {
            tiny_,
            _A1m1,
            _A2m1,
            _A3c,
            _A4,
            _B11,
            _B21,
            _B31,
            _B41,
            _C1a,
            _C1pa,
            _comg1,
            _C2a,
            _C3a,
            _C4a,
            _b,
            _c2,
            _calp0,
            _csig1,
            _ctau1,
            _dn1,
            _f1,
            _k2,
            _salp0,
            _somg1,
            _ssig1,
            _stau1,
            _a,
            _a13,
            azi1,
            calp1,
            caps,
            f,
            lat1,
            lon1,
            _s13,
            salp1,
        }
    }

    /// returns (a12, lat2, lon2, azi2, s12, m12, M12, M21, S12)
    pub fn _gen_position(
        &self,
        arcmode: bool,
        s12_a12: f64,
        outmask: u64,
    ) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64) {
        let mut a12 = std::f64::NAN;
        let mut lat2 = std::f64::NAN;
        let mut lon2 = std::f64::NAN;
        let mut azi2 = std::f64::NAN;
        let mut s12 = std::f64::NAN;
        let mut m12 = std::f64::NAN;
        let mut M12 = std::f64::NAN;
        let mut M21 = std::f64::NAN;
        let mut S12 = std::f64::NAN;
        let outmask = outmask & (self.caps & caps::OUT_MASK);
        if !(arcmode || (self.caps & (caps::OUT_MASK & caps::DISTANCE_IN) != 0)) {
            return (a12, lat2, lon2, azi2, s12, m12, M12, M21, S12);
        }

        let mut B12 = 0.0;
        let mut AB1 = 0.0;
        let mut sig12: f64;
        let mut ssig12: f64;
        let mut csig12: f64;
        let mut ssig2: f64;
        let mut csig2: f64;
        if arcmode {
            sig12 = s12_a12.to_radians();
            let res = geomath::sincosd(s12_a12);
            ssig12 = res.0;
            csig12 = res.1;
        } else {
            // tau12 = s12_a12 / (self._b * (1 + self._A1m1))
            let tau12 = s12_a12 / (self._b * (1.0 + self._A1m1));

            let s = tau12.sin();
            let c = tau12.cos();

            B12 = -geomath::sin_cos_series(
                true,
                self._stau1 * c + self._ctau1 * s,
                self._ctau1 * c - self._stau1 * s,
                &self._C1pa,
            );
            sig12 = tau12 - (B12 - self._B11);
            ssig12 = sig12.sin();
            csig12 = sig12.cos();
            if self.f.abs() > 0.01 {
                ssig2 = self._ssig1 * csig12 + self._csig1 * ssig12;
                csig2 = self._csig1 * csig12 - self._ssig1 * ssig12;
                B12 = geomath::sin_cos_series(true, ssig2, csig2, &self._C1a);
                let serr = (1.0 + self._A1m1) * (sig12 + (B12 - self._B11)) - s12_a12 / self._b;
                sig12 -= serr / (1.0 + self._k2 * geomath::sq(ssig2)).sqrt();
                ssig12 = sig12.sin();
                csig12 = sig12.cos();
            }
        };
        ssig2 = self._ssig1 * csig12 + self._csig1 * ssig12;
        csig2 = self._csig1 * csig12 - self._ssig1 * ssig12;
        let dn2 = (1.0 + self._k2 * geomath::sq(ssig2)).sqrt();
        if outmask & (caps::DISTANCE | caps::REDUCEDLENGTH | caps::GEODESICSCALE) != 0 {
            if arcmode || self.f.abs() > 0.01 {
                B12 = geomath::sin_cos_series(true, ssig2, csig2, &self._C1a);
            }
            AB1 = (1.0 + self._A1m1) * (B12 - self._B11);
        }

        let sbet2 = self._calp0 * ssig2;
        let mut cbet2 = self._salp0.hypot(self._calp0 * csig2);
        if cbet2 == 0.0 {
            cbet2 = self.tiny_;
            csig2 = self.tiny_;
        }
        let salp2 = self._salp0;
        let calp2 = self._calp0 * csig2;
        if outmask & caps::DISTANCE != 0 {
            s12 = if arcmode {
                self._b * ((1.0 + self._A1m1) * sig12 + AB1)
            } else {
                s12_a12
            }
        }
        if outmask & caps::LONGITUDE != 0 {
            let somg2 = self._salp0 * ssig2;
            let comg2 = csig2;
            let E = 1.0_f64.copysign(self._salp0);
            let omg12 = if outmask & caps::LONG_UNROLL != 0 {
                E * (sig12 - (ssig2.atan2(csig2) - self._ssig1.atan2(self._csig1))
                    + ((E * somg2).atan2(comg2) - (E * self._somg1).atan2(self._comg1)))
            } else {
                (somg2 * self._comg1 - comg2 * self._somg1)
                    .atan2(comg2 * self._comg1 + somg2 * self._somg1)
            };
            let lam12 = omg12
                + self._A3c
                    * (sig12
                        + (geomath::sin_cos_series(true, ssig2, csig2, &self._C3a) - self._B31));
            let lon12 = lam12.to_degrees();
            lon2 = if outmask & caps::LONG_UNROLL != 0 {
                self.lon1 + lon12
            } else {
                geomath::ang_normalize(
                    geomath::ang_normalize(self.lon1) + geomath::ang_normalize(lon12),
                )
            };
        };

        if outmask & caps::LATITUDE != 0 {
            lat2 = geomath::atan2d(sbet2, self._f1 * cbet2);
        }
        if outmask & caps::AZIMUTH != 0 {
            azi2 = geomath::atan2d(salp2, calp2);
        }
        if outmask & (caps::REDUCEDLENGTH | caps::GEODESICSCALE) != 0 {
            let B22 = geomath::sin_cos_series(true, ssig2, csig2, &self._C2a);
            let AB2 = (1.0 + self._A2m1) * (B22 - self._B21);
            let J12 = (self._A1m1 - self._A2m1) * sig12 + (AB1 - AB2);
            if outmask & caps::REDUCEDLENGTH != 0 {
                m12 = self._b
                    * ((dn2 * (self._csig1 * ssig2) - self._dn1 * (self._ssig1 * csig2))
                        - self._csig1 * csig2 * J12);
            }
            if outmask & caps::GEODESICSCALE != 0 {
                let t =
                    self._k2 * (ssig2 - self._ssig1) * (ssig2 + self._ssig1) / (self._dn1 + dn2);
                M12 = csig12 + (t * ssig2 - csig2 * J12) * self._ssig1 / self._dn1;
                M21 = csig12 - (t * self._ssig1 - self._csig1 * J12) * ssig2 / dn2;
            }
        }
        if outmask & caps::AREA != 0 {
            let B42 = geomath::sin_cos_series(false, ssig2, csig2, &self._C4a);
            let salp12: f64;
            let calp12: f64;
            if self._calp0 == 0.0 || self._salp0 == 0.0 {
                salp12 = salp2 * self.calp1 - calp2 * self.salp1;
                calp12 = calp2 * self.calp1 + salp2 * self.salp1;
            } else {
                salp12 = self._calp0
                    * self._salp0
                    * (if csig12 <= 0.0 {
                        self._csig1 * (1.0 - csig12) + ssig12 * self._ssig1
                    } else {
                        ssig12 * (self._csig1 * ssig12 / (1.0 + csig12) + self._ssig1)
                    });
                calp12 = geomath::sq(self._salp0) + geomath::sq(self._calp0) * self._csig1 * csig2;
            }
            S12 = self._c2 * salp12.atan2(calp12) + self._A4 * (B42 - self._B41);
        }
        a12 = if arcmode { s12_a12 } else { sig12.to_degrees() };
        (a12, lat2, lon2, azi2, s12, m12, M12, M21, S12)
    }

    // not currently used, but maybe some day
    #[allow(dead_code)]
    pub fn Position(&self, s12: f64, outmask: Option<u64>) -> HashMap<String, f64> {
        let outmask = match outmask {
            Some(outmask) => outmask,
            None => caps::STANDARD,
        };
        let mut result: HashMap<String, f64> = HashMap::new();
        result.insert("lat1".to_string(), self.lat1);
        result.insert("azi1".to_string(), self.azi1);
        result.insert("s12".to_string(), s12);
        let lon1 = if outmask & caps::LONG_UNROLL != 0 {
            self.lon1
        } else {
            geomath::ang_normalize(self.lon1)
        };
        result.insert("lon1".to_string(), lon1);

        let (a12, lat2, lon2, azi2, _s12, m12, M12, M21, S12) =
            self._gen_position(false, s12, outmask);
        let outmask = outmask & caps::OUT_MASK;
        result.insert("a12".to_string(), a12);
        if outmask & caps::LATITUDE != 0 {
            result.insert("lat2".to_string(), lat2);
        }
        if outmask & caps::LONGITUDE != 0 {
            result.insert("lon2".to_string(), lon2);
        }
        if outmask & caps::AZIMUTH != 0 {
            result.insert("azi2".to_string(), azi2);
        }
        if outmask & caps::REDUCEDLENGTH != 0 {
            result.insert("m12".to_string(), m12);
        }
        if outmask & caps::GEODESICSCALE != 0 {
            result.insert("M12".to_string(), M12);
            result.insert("M21".to_string(), M21);
        }
        if outmask & caps::AREA != 0 {
            result.insert("S12".to_string(), S12);
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use geodesic::Geodesic;

    #[test]
    fn test_gen_position() {
        let geod = Geodesic::wgs84();
        let gl = GeodesicLine::new(&geod, 0.0, 0.0, 10.0, None, None, None);
        let res = gl._gen_position(false, 150.0, 3979);
        assert_eq!(res.0, 0.0013520059461334633);
        assert_eq!(res.1, 0.0013359451088740494);
        assert_eq!(res.2, 0.00023398621812867812);
        assert_eq!(res.3, 10.000000002727887);
        assert_eq!(res.4, 150.0);
        assert!(res.5.is_nan());
        assert!(res.6.is_nan());
        assert!(res.7.is_nan());
        assert!(res.8.is_nan());
    }

    #[test]
    fn test_init() {
        let geod = Geodesic::wgs84();
        let gl = GeodesicLine::new(&geod, 0.0, 0.0, 0.0, None, None, None);
        assert_eq!(gl._a, 6378137.0);
        assert_eq!(gl.f, 0.0033528106647474805);
        assert_eq!(gl._b, 6356752.314245179);
        assert_eq!(gl._c2, 40589732499314.76);
        assert_eq!(gl._f1, 0.9966471893352525);
        assert_eq!(gl.caps, 36747);
        assert_eq!(gl.lat1, 0.0);
        assert_eq!(gl.lon1, 0.0);
        assert_eq!(gl.azi1, 0.0);
        assert_eq!(gl.salp1, 0.0);
        assert_eq!(gl.calp1, 1.0);
        assert_eq!(gl._dn1, 1.0);
        assert_eq!(gl._salp0, 0.0);
        assert_eq!(gl._calp0, 1.0);
        assert_eq!(gl._ssig1, 0.0);
        assert_eq!(gl._somg1, 0.0);
        assert_eq!(gl._csig1, 1.0);
        assert_eq!(gl._comg1, 1.0);
        assert_eq!(gl._k2, geod._ep2);
        assert!(gl._s13.is_nan());
        assert!(gl._a13.is_nan());
    }
}
