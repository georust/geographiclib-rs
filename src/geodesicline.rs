use crate::geodesic;
use crate::geodesiccapability as caps;
use crate::geomath;

#[derive(Debug)]
struct GeodesicLine {
    _A1m1: f64,
    _A2m1: f64,
    _A3c: f64,
    _A4: f64,
    _B11: f64,
    _B21: f64,
    _B31: f64,
    _B41: f64,
    _C1a: Vec<f64>,
    _C1pa: Vec<f64>,
    _C2a: Vec<f64>,
    _C3a: Vec<f64>,
    _C4a: Vec<f64>,
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
    a13: f64,
    a: f64,
    azi1: f64,
    calp1: f64,
    caps: u64,
    f: f64,
    lat1: f64,
    lon1: f64,
    s13: f64,
    salp1: f64,
}

impl GeodesicLine {
    pub fn new(
        geod: geodesic::Geodesic,
        lat1: f64,
        lon1: f64,
        azi1: f64,
        caps: u64,
        salp1: f64,
        calp1: f64,
    ) -> Self {
        let a = geod.a;
        let f = geod.f;
        let _b = geod._b;
        let _c2 = geod._c2;
        let _f1 = geod._f1;
        let caps = caps | caps::LATITUDE | caps::AZIMUTH | caps::LONG_UNROLL;
        let lat1 = geomath::lat_fix(lat1);
        let lon1 = lon1;
        let (azi1, salp1, calp1) = if salp1.is_nan() || calp1.is_nan() {
            let res = geomath::sincosd(geomath::ang_round(azi1));
            (geomath::ang_normalize(azi1), res.0, res.1)
        } else {
            (azi1, salp1, calp1)
        };

        let (mut sbet1, cbet1) = geomath::sincosd(geomath::ang_round(lat1));
        sbet1 *= _f1;
        let (sbet1, mut cbet1) = geomath::norm(sbet1, cbet1);
        cbet1 = geod.tiny_.max(cbet1);
        let _dn1 = (1.0 + geod._ep2 * geomath::sq(sbet1)).sqrt();
        let _salp0 = salp1 * cbet1;
        let _calp0 = calp1.hypot(salp1 * sbet1);
        let _ssig1 = sbet1;
        let _somg1 = _salp0 * sbet1;
        let _csig1 = if sbet1 != 0.0 || calp1 != 0.0 {
            cbet1 * calp1
        } else {
            1.0
        };
        let _comg1 = _csig1;
        let (_ssig1, _csig1) = geomath::norm(_ssig1, _csig1);
        let _k2 = geomath::sq(_calp0) * geod._ep2;
        let eps = _k2 / (2.0 * (1.0 + (1.0 + _k2).sqrt()) + _k2);

        let mut _A1m1 = 0.0;
        let mut _C1a = (0..geod.GEODESIC_ORDER).map(|x| x as f64).collect();
        let mut _B11 = 0.0;
        let mut _stau1 = 0.0;
        let mut _ctau1 = 0.0;
        if caps & caps::CAP_C1 != 0 {
            _A1m1 = geomath::_A1m1f(eps, geod.GEODESIC_ORDER);
            geomath::_C1f(eps, &mut _C1a, geod.GEODESIC_ORDER);
            _B11 = geomath::sin_cos_series(true, _ssig1, _csig1, _C1a.clone());
            let s = _B11.sin();
            let c = _B11.cos();
            _stau1 = _ssig1 * c + _csig1 * s;
            _ctau1 = _csig1 * c - _ssig1 * s;
        }

        let mut _C1pa = (0..=geod.GEODESIC_ORDER).map(|x| x as f64).collect();
        if caps & caps::CAP_C1p != 0 {
            geomath::_C1pf(eps, &mut _C1pa, geod.GEODESIC_ORDER);
        }

        let mut _A2m1 = 0.0;
        let mut _C2a = (0..=geod.GEODESIC_ORDER).map(|x| x as f64).collect();
        let mut _B21 = 0.0;
        if caps & caps::CAP_C2 != 0 {
            _A2m1 = geomath::_A2m1f(eps, geod.GEODESIC_ORDER);
            geomath::_C2f(eps, &mut _C2a, geod.GEODESIC_ORDER);
            _B21 = geomath::sin_cos_series(true, _ssig1, _csig1, _C2a.clone());
        }

        let mut _C3a = (0..geod.GEODESIC_ORDER).map(|x| x as f64).collect();
        let mut _A3c = 0.0;
        let mut _B31 = 0.0;
        if caps & caps::CAP_C3 != 0 {
            geod._C3f(eps, &mut _C3a);
            _A3c = -f * _salp0 * geod._A3f(eps);
            _B31 = geomath::sin_cos_series(true, _ssig1, _csig1, _C3a.clone());
        }

        let mut _C4a = (0..geod.GEODESIC_ORDER).map(|x| x as f64).collect();
        let mut _A4 = 0.0;
        let mut _B41 = 0.0;
        if caps & caps::CAP_C4 != 0 {
            geod._C4f(eps, &mut _C4a);
            _A4 = geomath::sq(a) * _calp0 * _salp0 * geod._e2;
            _B41 = geomath::sin_cos_series(false, _ssig1, _csig1, _C4a.clone());
        }

        let s13 = std::f64::NAN;
        let a13 = std::f64::NAN;

        GeodesicLine {
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
            a,
            a13,
            azi1,
            calp1,
            caps,
            f,
            lat1,
            lon1,
            s13,
            salp1,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use geodesic::{Geodesic, WGS84_A, WGS84_F};

    #[test]
    fn test_init() {
        let geod = Geodesic::new(WGS84_A, WGS84_F);
        let gl = GeodesicLine::new(geod, 0.0, 0.0, 0.0, 0, 0.0, 0.0);
        assert_eq!(gl.a, 6378137.0);
        assert_eq!(gl.f, 0.0033528106647474805);
        assert_eq!(gl._b, 6356752.314245179);
        assert_eq!(gl._c2, 40589732499314.76);
        assert_eq!(gl._f1, 0.9966471893352525);
        assert_eq!(gl.caps, 33408);
        assert_eq!(gl.lat1, 0.0);
        assert_eq!(gl.lon1, 0.0);
        assert_eq!(gl.azi1, 0.0);
        assert_eq!(gl.salp1, 0.0);
        assert_eq!(gl.calp1, 0.0);
        assert_eq!(gl._dn1, 1.0);
        assert_eq!(gl._salp0, 0.0);
        assert_eq!(gl._calp0, 0.0);
        assert_eq!(gl._ssig1, 0.0);
        assert_eq!(gl._somg1, 0.0);
        assert_eq!(gl._csig1, 1.0);
        assert_eq!(gl._comg1, 1.0);
        assert_eq!(gl._k2, 0.0);
        assert_eq!(gl.s13.is_nan(), true);
        assert_eq!(gl.a13.is_nan(), true);
    }
}
