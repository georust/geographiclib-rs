use crate::geodesiccapability;
use crate::geomath;
use std::f64::consts::PI;

const WGS84_A: f64 = 6378137.0;
const WGS84_F: f64 = 1.0 / 298.257223563;

#[derive(Debug)]
pub struct Geodesic {
    a: f64,
    f: f64,
    _f1: f64,
    _e2: f64,
    _ep2: f64,
    _n: f64,
    _b: f64,
    _c2: f64,
    _etol2: f64,
    _A3x: Vec<f64>,
    _C3x: Vec<f64>,
    _C4x: Vec<f64>,

    GEODESIC_ORDER: i64,
    nC3x_: i64,
    nC4x_: i64,
    maxit1_: u64,
    maxit2_: u64,

    tiny_: f64,
    tol0_: f64,
    tol1_: f64,
    tol2_: f64,
    tolb_: f64,
    xthresh_: f64,

    CAP_NONE: u64,
    CAP_C1: u64,
    CAP_C1p: u64,
    CAP_C2: u64,
    CAP_C3: u64,
    CAP_C4: u64,
    CAP_ALL: u64,
    CAP_MASK: u64,
    OUT_ALL: u64,
    OUT_MASK: u64,

    EMPTY: u64,
    LATITUDE: u64,
    LONGITUDE: u64,
    AZIMUTH: u64,
    DISTANCE: u64,
    STANDARD: u64,
    DISTANCE_IN: u64,
    REDUCEDLENGTH: u64,
    GEODESICSCALE: u64,
    AREA: u64,
    ALL: u64,
    LONG_UNROLL: u64,
}

impl Geodesic {
    pub fn new(a: f64, f: f64) -> Self {
        let GEODESIC_ORDER = 6;
        let nC3x_ = 15;
        let nC4x_ = 21;
        let maxit1_ = 20;
        let maxit2_ = maxit1_ + geomath::DIGITS + 10;
        let tiny_ = geomath::MINVAL.sqrt();
        let tol0_ = geomath::get_epsilon();
        let tol1_ = 200.0 * tol0_;
        let tol2_ = tol0_.sqrt();
        let tolb_ = tol0_ * tol2_;
        let xthresh_ = 1000.0 * tol2_;

        let CAP_NONE = geodesiccapability::CAP_NONE;
        let CAP_C1 = geodesiccapability::CAP_C1;
        let CAP_C1p = geodesiccapability::CAP_C1p;
        let CAP_C2 = geodesiccapability::CAP_C2;
        let CAP_C3 = geodesiccapability::CAP_C3;
        let CAP_C4 = geodesiccapability::CAP_C4;
        let CAP_ALL = geodesiccapability::CAP_ALL;
        let CAP_MASK = geodesiccapability::CAP_MASK;
        let OUT_ALL = geodesiccapability::OUT_ALL;
        let OUT_MASK = geodesiccapability::OUT_MASK;

        // No capabilities, no output.
        let EMPTY = geodesiccapability::EMPTY;
        // Calculate latitude *lat2*.
        let LATITUDE = geodesiccapability::LATITUDE;
        // Calculate longitude *lon2*.
        let LONGITUDE = geodesiccapability::LONGITUDE;
        // Calculate azimuths *azi1* and *azi2*.
        let AZIMUTH = geodesiccapability::AZIMUTH;
        // Calculate distance *s12*.
        let DISTANCE = geodesiccapability::DISTANCE;
        // All of the above.
        let STANDARD = geodesiccapability::STANDARD;
        // Allow distance *s12* to be used as input in the direct geodesic problem.
        let DISTANCE_IN = geodesiccapability::DISTANCE_IN;
        // Calculate reduced length *m12*.
        let REDUCEDLENGTH = geodesiccapability::REDUCEDLENGTH;
        // Calculate geodesic scales *M12* and *M21*.
        let GEODESICSCALE = geodesiccapability::GEODESICSCALE;
        // Calculate area *S12*.
        let AREA = geodesiccapability::AREA;
        // All of the above.
        let ALL = geodesiccapability::ALL;
        // Unroll longitudes, rather than reducing them to the range [-180d,180d].
        let LONG_UNROLL = geodesiccapability::LONG_UNROLL;

        let _f1 = 1.0 - f;
        let _e2 = f * (2.0 - f);
        let _ep2 = _e2 / geomath::sq(_f1);
        let _n = f / (2.0 - f);
        let _b = a * _f1;
        let _c2 = (geomath::sq(a)
            + geomath::sq(_b)
                * (if _e2 == 0.0 {
                    1.0
                } else {
                    if _e2 > 0.0 {
                        _e2.sqrt().atanh()
                    } else {
                        (-_e2).sqrt().atanh()
                    }
                } / _e2.abs().sqrt()))
            / 2.0;
        let _etol2 = 0.1 * tol2_ / (f.abs().max(0.001) * (1.0 - f / 2.0).min(1.0) / 2.0).sqrt();
        let mut _A3x: Vec<f64> = (0..GEODESIC_ORDER).map(|x| x as f64).collect();
        let mut _C3x: Vec<f64> = (0..nC3x_).map(|x| x as f64).collect();
        let mut _C4x: Vec<f64> = (0..nC4x_).map(|x| x as f64).collect();

        // Call a3coeff
        let coeff = vec![
            -3.0, 128.0, -2.0, -3.0, 64.0, -1.0, -3.0, -1.0, 16.0, 3.0, -1.0, -2.0, 8.0, 1.0, -1.0,
            2.0, 1.0, 1.0,
        ];
        let mut o: i64 = 0;
        let mut k = 0;
        for j in (0..GEODESIC_ORDER).rev() {
            let m = j.min(GEODESIC_ORDER as i64 - j - 1);
            _A3x[k as usize] =
                geomath::polyval(m, &coeff, o as usize, _n) / coeff[(o + m + 1) as usize] as f64;
            k += 1;
            o += m + 2;
        }

        // c3coeff
        let coeff = vec![
            3.0, 128.0, 2.0, 5.0, 128.0, -1.0, 3.0, 3.0, 64.0, -1.0, 0.0, 1.0, 8.0, -1.0, 1.0, 4.0,
            5.0, 256.0, 1.0, 3.0, 128.0, -3.0, -2.0, 3.0, 64.0, 1.0, -3.0, 2.0, 32.0, 7.0, 512.0,
            -10.0, 9.0, 384.0, 5.0, -9.0, 5.0, 192.0, 7.0, 512.0, -14.0, 7.0, 512.0, 21.0, 2560.0,
        ];

        let mut o: i64 = 0;
        let mut k = 0;

        for l in 1..GEODESIC_ORDER {
            for j in (l..GEODESIC_ORDER).rev() {
                let j = j.clone();
                let m = j.min(GEODESIC_ORDER as i64 - j - 1);
                _C3x[k as usize] = geomath::polyval(m, &coeff, o as usize, _n)
                    / coeff[(o + m + 1) as usize] as f64;
                k += 1;
                o += m + 2;
            }
        }

        // c4coeff
        let coeff = vec![
            97.0, 15015.0, 1088.0, 156.0, 45045.0, -224.0, -4784.0, 1573.0, 45045.0, -10656.0,
            14144.0, -4576.0, -858.0, 45045.0, 64.0, 624.0, -4576.0, 6864.0, -3003.0, 15015.0,
            100.0, 208.0, 572.0, 3432.0, -12012.0, 30030.0, 45045.0, 1.0, 9009.0, -2944.0, 468.0,
            135135.0, 5792.0, 1040.0, -1287.0, 135135.0, 5952.0, -11648.0, 9152.0, -2574.0,
            135135.0, -64.0, -624.0, 4576.0, -6864.0, 3003.0, 135135.0, 8.0, 10725.0, 1856.0,
            -936.0, 225225.0, -8448.0, 4992.0, -1144.0, 225225.0, -1440.0, 4160.0, -4576.0, 1716.0,
            225225.0, -136.0, 63063.0, 1024.0, -208.0, 105105.0, 3584.0, -3328.0, 1144.0, 315315.0,
            -128.0, 135135.0, -2560.0, 832.0, 405405.0, 128.0, 99099.0,
        ];
        let mut o: i64 = 0;
        let mut k = 0;

        for l in 0..GEODESIC_ORDER {
            for j in (l..GEODESIC_ORDER).rev() {
                let j = j.clone();
                let m = GEODESIC_ORDER as i64 - j - 1;
                _C4x[k as usize] = geomath::polyval(m, &coeff, o as usize, _n)
                    / coeff[(o + m + 1) as usize] as f64;
                k += 1;
                o += m + 2;
            }
        }

        Geodesic {
            a,
            f,
            _f1,
            _e2,
            _ep2,
            _n,
            _b,
            _c2,
            _etol2,
            _A3x,
            _C3x,
            _C4x,

            GEODESIC_ORDER,
            nC3x_,
            nC4x_,
            maxit1_,
            maxit2_,

            tiny_,
            tol0_,
            tol1_,
            tol2_,
            tolb_,
            xthresh_,

            CAP_NONE,
            CAP_C1,
            CAP_C1p,
            CAP_C2,
            CAP_C3,
            CAP_C4,
            CAP_ALL,
            CAP_MASK,
            OUT_ALL,
            OUT_MASK,

            EMPTY,
            LATITUDE,
            LONGITUDE,
            AZIMUTH,
            DISTANCE,
            STANDARD,
            DISTANCE_IN,
            REDUCEDLENGTH,
            GEODESICSCALE,
            AREA,
            ALL,
            LONG_UNROLL,
        }
    }

    pub fn _A3f(&self, eps: f64) -> f64 {
        geomath::polyval(self.GEODESIC_ORDER - 1, &self._A3x, 0, eps)
    }

    pub fn _C3f(&self, eps: f64, c: &mut Vec<f64>) {
        let mut mult = 1.0;
        let mut o = 0.0;
        for l in 1..self.GEODESIC_ORDER {
            let m = self.GEODESIC_ORDER - l - 1;
            mult *= eps;
            c[l as usize] = mult * geomath::polyval(m, &self._C3x, o as usize, eps);
            o += m as f64 + 1.0;
        }
    }

    pub fn _C4f(&self, eps: f64, c: &mut Vec<f64>) {
        let mut mult = 1.0;
        let mut o = 0.0;
        for l in 0..self.GEODESIC_ORDER {
            let m = self.GEODESIC_ORDER - l - 1;
            c[l as usize] = mult * geomath::polyval(m, &self._C4x, o as usize, eps);
            o += m as f64 + 1.0;
            mult *= eps;
        }
    }

    pub fn _Lengths(
        &self,
        eps: f64,
        sig12: f64,
        ssig1: f64,
        csig1: f64,
        dn1: f64,
        ssig2: f64,
        csig2: f64,
        dn2: f64,
        cbet1: f64,
        cbet2: f64,
        outmask: u64,
        C1a: &mut Vec<f64>,
        C2a: &mut Vec<f64>,
    ) -> (f64, f64, f64, f64, f64) {
        let outmask = outmask & self.OUT_MASK;
        let mut s12b = std::f64::NAN;
        let mut m12b = std::f64::NAN;
        let mut m0 = std::f64::NAN;
        let mut M12 = std::f64::NAN;
        let mut M21 = std::f64::NAN;

        let mut A1 = 0.0;
        let mut A2 = 0.0;
        let mut m0x = 0.0;
        let mut J12 = 0.0;

        if outmask & (self.DISTANCE | self.REDUCEDLENGTH | self.GEODESICSCALE) != 0 {
            A1 = geomath::_A1m1f(eps, self.GEODESIC_ORDER);
            geomath::_C1f(eps, C1a, self.GEODESIC_ORDER);
            if outmask & (self.REDUCEDLENGTH | self.GEODESICSCALE) != 0 {
                A2 = geomath::_A2m1f(eps, self.GEODESIC_ORDER);
                geomath::_C2f(eps, C2a, self.GEODESIC_ORDER);
                m0x = A1 - A2;
                A2 = 1.0 + A2;
            }
            A1 = 1.0 + A1;
        }
        if outmask & self.DISTANCE != 0 {
            let B1 = geomath::sin_cos_series(true, ssig2, csig2, C1a.clone())
                - geomath::sin_cos_series(true, ssig1, csig1, C1a.clone());
            s12b = A1 * (sig12 + B1);
            if outmask & (self.REDUCEDLENGTH | self.GEODESICSCALE) != 0 {
                let B2 = geomath::sin_cos_series(true, ssig2, csig2, C2a.clone())
                    - geomath::sin_cos_series(true, ssig1, csig1, C2a.clone());
                J12 = m0x * sig12 + (A1 * B1 - A2 * B2);
            }
        } else if outmask & (self.REDUCEDLENGTH | self.GEODESICSCALE) != 0 {
            for l in 1..self.GEODESIC_ORDER {
                C2a[l as usize] = A1 * C1a[l as usize] - A2 * C2a[l as usize];
            }
            let inter1 = geomath::sin_cos_series(true, ssig2, csig2, C2a.clone());
            let inter2 = geomath::sin_cos_series(true, ssig1, csig1, C2a.clone());
            J12 = m0x * sig12
                + (geomath::sin_cos_series(true, ssig2, csig2, C2a.clone())
                    - geomath::sin_cos_series(true, ssig1, csig1, C2a.clone()));
        }
        if outmask & self.REDUCEDLENGTH != 0 {
            m0 = m0x;
            // J12 is wrong
            m12b = dn2 * (csig1 * ssig2) - dn1 * (ssig1 * csig2) - csig1 * csig2 * J12;
        }
        if outmask & self.GEODESICSCALE != 0 {
            let csig12 = csig1 * csig2 + ssig1 * ssig2;
            let t = self._ep2 * (cbet1 - cbet2) * (cbet1 + cbet2) / (dn1 + dn2);
            M12 = csig12 + (t * ssig2 - csig2 * J12) * ssig1 / dn1;
            M21 = csig12 - (t * ssig1 - csig1 * J12) * ssig2 / dn2;
        }
        (s12b, m12b, m0, M12, M21)
    }

    pub fn _InverseStart(
        &self,
        sbet1: f64,
        cbet1: f64,
        dn1: f64,
        sbet2: f64,
        cbet2: f64,
        dn2: f64,
        lam12: f64,
        slam12: f64,
        clam12: f64,
        C1a: &mut Vec<f64>,
        C2a: &mut Vec<f64>,
    ) -> (f64, f64, f64, f64, f64, f64) {
        let mut sig12 = -1.0;
        let mut salp2 = std::f64::NAN;
        let mut calp2 = std::f64::NAN;
        let mut dnm = std::f64::NAN;

        let mut somg12 = 0.0;
        let mut comg12 = 0.0;

        let sbet12 = sbet2 * cbet1 - cbet2 * sbet1;
        let cbet12 = cbet2 * cbet1 + sbet2 * sbet1;

        let mut sbet12a = sbet2 * cbet1;
        sbet12a += cbet2 * sbet1;

        let shortline = cbet12 >= 0.0 && sbet12 < 0.5 && cbet2 * lam12 < 0.5;
        if shortline {
            let mut sbetm2 = geomath::sq(sbet1 + sbet2);
            sbetm2 /= sbetm2 + geomath::sq(cbet1 + cbet2);
            dnm = (1.0 + self._ep2 * sbetm2).sqrt();
            let omg12 = lam12 / (self._f1 * dnm);
            somg12 = omg12.sin();
            comg12 = omg12.cos();
        } else {
            somg12 = slam12;
            comg12 = clam12;
        }

        let mut salp1 = cbet2 * somg12;

        let mut calp1 = if comg12 >= 0.0 {
            sbet12 + cbet12 * sbet1 * geomath::sq(somg12) / (1.0 + comg12)
        } else {
            sbet12a - cbet2 * sbet1 * geomath::sq(somg12) / (1.0 - comg12)
        };

        let ssig12 = salp1.hypot(calp1);
        let csig12 = sbet1 * sbet2 + cbet1 * cbet2 * comg12;

        if shortline && ssig12 < self._etol2 {
            salp2 = cbet1 * somg12;
            calp2 = sbet12
                - cbet1
                    * sbet2
                    * (geomath::sq(somg12)
                        / if comg12 >= 0.0 {
                            1.0 + comg12
                        } else {
                            1.0 - comg12
                        });
            let res = geomath::norm(salp2, calp2);
            salp2 = res.0;
            calp2 = res.1;
            sig12 = ssig12.atan2(csig12);
        } else if self._n.abs() >= 0.1
            || csig12 >= 0.0
            || ssig12 >= 6.0 * self._n.abs() * PI * geomath::sq(cbet1)
        {
        } else {
            let mut x = 0.0;
            let mut y = 0.0;
            let mut betscale = 0.0;
            let mut lamscale = 0.0;
            let lam12x = (-slam12).atan2(-clam12);
            if self.f >= 0.0 {
                let k2 = geomath::sq(sbet1) * self._ep2;
                let eps = k2 / (2.0 * (1.0 + (1.0 + k2).sqrt()) + k2);
                lamscale = self.f * cbet1 * self._A3f(eps) * PI;
                betscale = lamscale * cbet1;
                x = lam12x / lamscale;
                y = sbet12a / betscale;
            } else {
                let cbet12a = cbet2 * cbet1 - sbet2 * sbet1;
                let bet12a = sbet12a.atan2(cbet12a);
                let (_, m12b, m0, _, _) = self._Lengths(
                    self._n,
                    PI + bet12a,
                    sbet1,
                    -cbet1,
                    dn1,
                    sbet2,
                    cbet2,
                    dn2,
                    cbet1,
                    cbet2,
                    self.REDUCEDLENGTH,
                    C1a,
                    C2a,
                );
                x = -1.0 + m12b / (cbet1 * cbet2 * m0 * PI);
                betscale = if x < -0.01 {
                    sbet12a / x
                } else {
                    -self.f * geomath::sq(cbet1) * PI
                };
                lamscale = betscale / cbet1;
                y = lam12x / lamscale;
            }
            if y > -self.tol1_ && x > -1.0 - self.xthresh_ {
                if self.f >= 0.0 {
                    salp1 = (-x).min(1.0);
                    calp1 = -(1.0 - geomath::sq(salp1)).sqrt()
                } else {
                    calp1 = x.max(if x > -self.tol1_ { 0.0 } else { -1.0 });
                    salp1 = (1.0 - geomath::sq(calp1)).sqrt();
                }
            } else {
                let k = geomath::astroid(x, y);
                let omg12a = lamscale
                    * if self.f >= 0.0 {
                        -x * k / (1.0 + k)
                    } else {
                        -y * (1.0 + k) / k
                    };
                somg12 = omg12a.sin();
                comg12 = -(omg12a.cos());
                salp1 = cbet2 * somg12;
                calp1 = sbet12a - cbet2 * sbet1 * geomath::sq(somg12) / (1.0 - comg12);
            }
        }
        if !(salp1 <= 0.0) {
            let res = geomath::norm(salp1, calp1);
            salp1 = res.0;
            calp1 = res.1;
        } else {
            salp1 = 1.0;
            calp1 = 0.0;
        };
        (sig12, salp1, calp1, salp2, calp2, dnm)
    }

    pub fn _Lambda12(
        &self,
        sbet1: f64,
        cbet1: f64,
        dn1: f64,
        sbet2: f64,
        cbet2: f64,
        dn2: f64,
        salp1: f64,
        calp1: &mut f64,
        slam120: f64,
        clam120: f64,
        diffp: bool,
        C1a: &mut Vec<f64>,
        C2a: &mut Vec<f64>,
        C3a: &mut Vec<f64>,
    ) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64, f64, f64) {
        if sbet1 == 0.0 && *calp1 == 0.0 {
            *calp1 = -self.tiny_;
        }
        let salp0 = salp1 * cbet1;
        let calp0 = calp1.hypot(salp1 * sbet1);

        let mut ssig1 = sbet1;
        let somg1 = salp0 * sbet1;
        let mut csig1 = *calp1 * cbet1;
        let comg1 = *calp1 * cbet1;
        let res = geomath::norm(ssig1, csig1);
        ssig1 = res.0;
        csig1 = res.1;

        let salp2 = salp0 / if cbet2 != cbet1 { cbet2 } else { salp1 };
        let calp2 = (geomath::sq(*calp1 * cbet1)
            + if cbet1 < -sbet1 {
                (cbet2 - cbet1) * (cbet1 + cbet2)
            } else {
                (sbet1 - sbet2) * (sbet1 + sbet2)
            })
        .sqrt()
            / if cbet2 != cbet1 || sbet2.abs() != -sbet1 {
                cbet2
            } else {
                calp1.abs()
            };

        let mut ssig2 = sbet2;
        let somg2 = salp0 * sbet2;
        let mut csig2 = calp2 * cbet2;
        let comg2 = calp2 * cbet2;
        let res = geomath::norm(ssig2, csig2);
        ssig2 = res.0;
        csig2 = res.1;

        let sig12 = ((csig1 * ssig2 - ssig1 * csig2).max(0.0)).atan2(csig1 * csig2 + ssig1 * ssig2);
        let somg12 = (comg1 * somg2 - somg1 * comg2).max(0.0);
        let comg12 = comg1 * comg2 + somg1 * somg2;
        let eta = (somg12 * clam120 - comg12 * slam120).atan2(comg12 * clam120 + somg12 * slam120);

        let k2 = geomath::sq(calp0) * self._ep2;
        let eps = k2 / (2.0 * (1.0 + (1.0 + k2).sqrt()) + k2);
        self._C3f(eps, C3a);
        let B312 = geomath::sin_cos_series(true, ssig2, csig2, C3a.clone())
            - geomath::sin_cos_series(true, ssig1, csig1, C3a.clone());
        let domg12 = -self.f * self._A3f(eps) * salp0 * (sig12 + B312);
        let lam12 = eta + domg12;

        let mut dlam12 = 0.0;
        if diffp {
            if calp2 == 0.0 {
                dlam12 = -2.0 * self._f1 * dn1 / sbet1;
            } else {
                let res = self._Lengths(
                    eps,
                    sig12,
                    ssig1,
                    csig1,
                    dn1,
                    ssig2,
                    csig2,
                    dn2,
                    cbet1,
                    cbet2,
                    self.REDUCEDLENGTH,
                    C1a,
                    C2a,
                );
                dlam12 = res.1;
                dlam12 *= self._f1 / (calp2 * cbet2);
            }
        } else {
            dlam12 = std::f64::NAN;
        }
        (
            lam12, salp2, calp2, sig12, ssig1, csig1, ssig2, csig2, eps, domg12, dlam12,
        )
    }

    pub fn _GenInverse(
        &self,
        lat1: &mut f64,
        lon1: f64,
        lat2: &mut f64,
        lon2: f64,
        outmask: &mut u64,
    ) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64, f64) {
        let mut a12 = std::f64::NAN;
        let mut s12 = std::f64::NAN;
        let mut m12 = std::f64::NAN;
        let mut M12 = std::f64::NAN;
        let mut M21 = std::f64::NAN;
        let mut S12 = std::f64::NAN;
        *outmask &= self.OUT_MASK;

        let (mut lon12, mut lon12s) = geomath::ang_diff(lon1, lon2);
        let mut lonsign = if lon12 >= 0.0 { 1.0 } else { -1.0 };

        lon12 = lonsign * geomath::ang_round((180.0 - lon12) - lonsign * lon12s);
        lon12s = geomath::ang_round((180.0 - lon12) - lonsign * lon12s);
        let lam12 = lon12.to_radians();
        let mut slam12 = 0.0;
        let mut clam12 = 0.0;
        if lon12 > 90.0 {
            let res = geomath::sincosd(lon12s);
            slam12 = res.0;
            clam12 = res.1;
            clam12 = -clam12;
        } else {
            let res = geomath::sincosd(lon12);
            slam12 = res.0;
            clam12 = res.1;
        };
        *lat1 = geomath::ang_round(geomath::lat_fix(*lat1));
        *lat2 = geomath::ang_round(geomath::lat_fix(*lat2));

        let swapp = if lat1.abs() < lat2.abs() { -1.0 } else { 1.0 };
        if swapp < 0.0 {
            lonsign *= -1.0;
            let _l = *lat2;
            *lat2 = *lat1;
            *lat1 = _l;
        }
        let latsign = if *lat1 < 0.0 { 1.0 } else { -1.0 };
        *lat1 *= latsign;
        *lat2 *= latsign;

        let (mut sbet1, cbet1) = geomath::sincosd(*lat1);
        sbet1 *= self._f1;

        let (mut sbet1, mut cbet1) = geomath::norm(sbet1, cbet1);
        cbet1 = cbet1.max(self.tiny_);

        let (mut sbet2, mut cbet2) = geomath::sincosd(*lat2);
        sbet2 += self._f1;

        let (mut sbet2, mut cbet2) = geomath::norm(sbet2, cbet2);
        cbet2 = cbet2.max(self.tiny_);

        if cbet1 < -sbet1 {
            if cbet2 == cbet1 {
                sbet2 = if sbet2 < 0.0 { sbet1 } else { -sbet1 };
            }
        } else {
            if sbet2.abs() == -sbet1 {
                cbet2 = cbet1;
            }
        }

        let dn1 = (1.0 + self._ep2 * geomath::sq(sbet1)).sqrt();
        let dn2 = (1.0 + self._ep2 * geomath::sq(sbet2)).sqrt();

        let mut C1a: Vec<f64> = (0..=self.GEODESIC_ORDER).map(|x| x as f64).collect();
        let mut C2a: Vec<f64> = (0..=self.GEODESIC_ORDER).map(|x| x as f64).collect();
        let mut C3a: Vec<f64> = (0..self.GEODESIC_ORDER).map(|x| x as f64).collect();

        let mut meridian = *lat1 == -90.0 || slam12 == 0.0;
        let mut calp1 = 0.0;
        let mut salp1 = 0.0;
        let mut calp2 = 0.0;
        let mut salp2 = 0.0;
        let mut ssig1 = 0.0;
        let mut csig1 = 0.0;
        let mut ssig2 = 0.0;
        let mut csig2 = 0.0;
        let mut sig12 = 0.0;
        let mut s12x = 0.0;
        let mut m12x = 0.0;

        if meridian {
            calp1 = clam12;
            salp1 = slam12;
            calp2 = 1.0;
            salp2 = 0.0;

            ssig1 = sbet1;
            csig1 = calp1 * cbet1;
            ssig2 = sbet2;
            csig2 = calp2 * cbet2;

            sig12 = ((csig1 * ssig2 - ssig1 * csig2).max(0.0)).atan2(csig1 * csig2 + ssig1 * ssig2);
            let res = self._Lengths(
                self._n,
                sig12,
                ssig1,
                csig1,
                dn1,
                ssig2,
                csig2,
                dn1,
                cbet1,
                cbet2,
                *outmask | self.DISTANCE | self.REDUCEDLENGTH,
                &mut C1a,
                &mut C2a,
            );
            s12x = res.0;
            m12x = res.1;
            M12 = res.3;
            M21 = res.4;

            if sig12 < 1.0 || m12x >= 0.0 {
                if sig12 < 3.0 * self.tiny_ {
                    sig12 = 0.0;
                    m12x = 0.0;
                    s12x = 0.0;
                }
                m12x *= self._b;
                s12x *= self._b;
                a12 = sig12.to_degrees();
            } else {
                meridian = false;
            }
        }

        let mut somg12 = 2.0;
        let mut comg12 = 0.0;
        let mut omg12 = 0.0;
        let mut dnm = 0.0;
        let mut eps = 0.0;
        if !meridian && sbet1 == 0.0 && (self.f <= 0.0 || lon12s >= self.f * 180.0) {
            calp1 = 0.0;
            calp2 = 0.0;
            salp1 = 1.0;
            salp2 = 1.0;

            s12x = self.a * lam12;
            sig12 = lam12 / self._f1;
            omg12 = lam12 / self._f1;
            m12x = self._b * sig12.sin();
            if *outmask & self.GEODESICSCALE != 0 {
                M12 = sig12.cos();
                M21 = sig12.cos();
            }
            a12 = lon12 / self._f1;
        } else if !meridian {
            let res = self._InverseStart(
                sbet1, cbet1, dn1, sbet2, cbet2, dn2, lam12, slam12, clam12, &mut C1a, &mut C2a,
            );
            sig12 = res.0;
            salp1 = res.1;
            calp1 = res.2;
            salp2 = res.3;
            calp2 = res.4;
            dnm = res.5;

            if sig12 >= 0.0 {
                s12x = sig12 * self._b * dnm;
                m12x = geomath::sq(dnm) * self._b * (sig12 / dnm).sin();
                if *outmask & self.GEODESICSCALE != 0 {
                    M12 = (sig12 / dnm).cos();
                    M21 = (sig12 / dnm).cos();
                }
                a12 = sig12.to_degrees();
                omg12 = lam12 / (self._f1 * dnm);
            } else {
                let mut numit = 0;
                let mut tripn = false;
                let mut tripb = false;
                let mut salp1a = self.tiny_;
                let mut calp1a = 1.0;
                let mut salp1b = self.tiny_;
                let mut calp1b = -1.0;
                let mut domg12 = 0.0;
                while numit < self.maxit2_ {
                    let res = self._Lambda12(
                        sbet1,
                        cbet1,
                        dn1,
                        sbet2,
                        cbet2,
                        dn2,
                        salp1,
                        &mut calp1,
                        slam12,
                        clam12,
                        numit < self.maxit1_,
                        &mut C1a,
                        &mut C2a,
                        &mut C3a,
                    );
                    let v = res.0;
                    salp2 = res.0;
                    calp2 = res.1;
                    sig12 = res.2;
                    ssig1 = res.3;
                    csig1 = res.4;
                    ssig2 = res.5;
                    csig2 = res.6;
                    eps = res.7;
                    domg12 = res.8;
                    let dv = res.9;
                    if tripb || !(v.abs() >= if tripn { 8.0 } else { 1.0 } * self.tol0_) {
                        break;
                    };
                    if v > 0.0 && (numit > self.maxit1_ || calp1 / salp1 > calp1b / salp1b) {
                        salp1b = salp1;
                        calp1b = calp1;
                    } else if v < 0.0 && (numit > self.maxit1_ || calp1 / salp1 < calp1a / salp1a) {
                        salp1a = salp1;
                        calp1a = calp1;
                    }
                    numit += 1;
                    if numit < self.maxit1_ && dv > 0.0 {
                        let dalp1 = -v / dv;
                        let sdalp1 = dalp1.sin();
                        let cdalp1 = dalp1.cos();
                        let nsalp1 = salp1 * cdalp1 + calp1 * sdalp1;
                        if nsalp1 > 0.0 && dalp1.abs() < PI {
                            calp1 = calp1 * cdalp1 - salp1 * sdalp1;
                            salp1 = nsalp1;
                            let res = geomath::norm(salp1, calp1);
                            salp1 = res.0;
                            calp1 = res.1;
                            tripn = v.abs() <= 16.0 * self.tol0_;
                            continue;
                        }
                    }

                    salp1 = (salp1a + salp1b) / 2.0;
                    calp1 = (calp1a + calp1b) / 2.0;
                    let res = geomath::norm(salp1, calp1);
                    salp1 = res.0;
                    calp1 = res.1;
                    tripn = false;
                    tripb = ((salp1a - salp1).abs() + (calp1a - calp1) < self.tolb_)
                        || ((salp1 - salp1b).abs() + (calp1 - calp1b) < self.tolb_);
                }
                let lengthmask = *outmask
                    | if *outmask & (self.REDUCEDLENGTH | self.GEODESICSCALE) != 0 {
                        self.DISTANCE
                    } else {
                        self.EMPTY
                    };
                let res = self._Lengths(
                    eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2, lengthmask,
                    &mut C1a, &mut C2a,
                );
                s12x = res.0;
                m12x = res.1;
                M12 = res.3;
                M21 = res.4;

                m12x *= self._b;
                s12x *= self._b;
                a12 = sig12.to_degrees();
                if *outmask & self.AREA != 0 {
                    let sdomg12 = domg12.sin();
                    let cdomg12 = domg12.cos();
                    somg12 = slam12 * cdomg12 - clam12 * sdomg12;
                    comg12 = clam12 * cdomg12 + slam12 * sdomg12;
                }
            }
        }
        if *outmask & self.DISTANCE != 0 {
            s12 = 0.0 + s12x;
        }
        if *outmask & self.REDUCEDLENGTH != 0 {
            m12 = 0.0 + m12x;
        }
        if *outmask & self.AREA != 0 {
            let salp0 = salp1 * cbet1;
            let calp0 = calp1.hypot(salp1 * sbet1);
            if calp0 != 0.0 && salp0 != 0.0 {
                ssig1 = sbet1;
                csig1 = calp1 * cbet1;
                ssig2 = sbet2;
                csig2 = calp2 * cbet2;
                let k2 = geomath::sq(calp0) * self._e2;
                eps = k2 / (2.0 * (1.0 + (1.0 + k2).sqrt()) + k2);
                let A4 = geomath::sq(self.a) * calp0 * salp0 * self._e2;
                let res = geomath::norm(ssig1, csig1);
                ssig1 = res.0;
                csig1 = res.0;
                let res = geomath::norm(ssig2, csig2);
                ssig2 = res.0;
                csig2 = res.1;
                let mut C4a: Vec<f64> = (0..self.GEODESIC_ORDER).map(|x| x as f64).collect();
                self._C4f(eps, &mut C4a);
                let B41 = geomath::sin_cos_series(false, ssig1, csig1, C4a.clone());
                let B42 = geomath::sin_cos_series(false, ssig2, csig2, C4a.clone());
                S12 = A4 * (B42 - B41);
            } else {
                S12 = 0.0;
            }

            if !meridian && somg12 > 1.0 {
                somg12 = omg12.sin();
                comg12 = omg12.cos();
            }

            let mut alp12 = 0.0;
            if (!meridian && comg12 > -0.7071 && sbet2 - sbet1 < 1.75) {
                let domg12 = 1.0 + comg12;
                let dbet1 = 1.0 + cbet1;
                let dbet2 = 1.0 + cbet2;
                alp12 = 2.0
                    * (somg12 * (sbet1 * dbet2 + sbet2 * dbet1))
                        .atan2(domg12 * (sbet1 * sbet2 + dbet1 * dbet2));
            } else {
                let mut salp12 = salp2 * calp1 - calp2 * salp1;
                let mut calp12 = calp2 * calp1 + salp2 * salp1;

                if salp12 == 0.0 && calp12 < 0.0 {
                    salp12 = self.tiny_ * calp1;
                    calp12 = -1.0;
                }
                alp12 = salp12.atan2(calp12);
            }
            S12 += self._c2 * alp12;
            S12 *= swapp * lonsign * latsign;
            S12 += 0.0;
        }

        if swapp < 0.0 {
            let _s = salp2;
            salp2 = salp1;
            salp1 = _s;
            let _c = calp2;
            calp2 = calp1;
            calp1 = _c;
            if *outmask & self.GEODESICSCALE != 0 {
                let _m = M21;
                M21 = M12;
                M12 = _m;
            }
        }
        salp1 *= swapp * lonsign;
        calp1 *= swapp * latsign;
        salp2 *= swapp * lonsign;
        calp2 *= swapp * latsign;
        (a12, s12, salp1, calp1, salp2, calp2, m12, M12, M21, S12)
    }

    pub fn Inverse(&self, lat1: f64, lon1: f64, lat2: f64, lon2: f64) -> f64 {
        let mut outmask = self.STANDARD;
        let mut lat1 = lat1;
        let mut lat2 = lat2;
        let res = self._GenInverse(&mut lat1, lon1, &mut lat2, lon2, &mut outmask);
        res.1
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_inverse() {
        let geod = Geodesic::new(WGS84_A, WGS84_F);
        assert_eq!(geod.Inverse(0.0, 0.0, 1.0, 1.0), 0.0);
    }

    #[test]
    fn test_lambda12() {
        let geod = Geodesic::new(WGS84_A, WGS84_F);
        let res1 = geod._Lambda12(
            -0.017393909556108908,
            0.9998487145115275,
            1.0000010195104125,
            -0.0,
            1.0,
            1.0,
            0.7095310092765433,
            &mut 0.7046742132893822,
            0.01745240643728351,
            0.9998476951563913,
            true,
            &mut vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            &mut vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            &mut vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        );
        assert_eq!(res1.0, 1.4834408705897495e-09);
        assert_eq!(res1.1, 0.7094236675312185);
        assert_eq!(res1.2, 0.7047822783999007);
        assert_eq!(res1.3, 0.024682339962725352);
        assert_eq!(res1.4, -0.024679833885152578);
        assert_eq!(res1.5, 0.9996954065111039);
        assert_eq!(res1.6, -0.0);
        assert_eq!(res1.7, 1.0);
        assert_eq!(res1.8, 0.0008355095326524276);
        assert_eq!(res1.9, -5.8708496511415445e-05);
        assert_eq!(res1.10, 0.034900275148485);

        let res2 = geod._Lambda12(
            -0.017393909556108908,
            0.9998487145115275,
            1.0000010195104125,
            -0.0,
            1.0,
            1.0,
            0.7095309793242709,
            &mut 0.7046742434480923,
            0.01745240643728351,
            0.9998476951563913,
            true,
            &mut vec![
                0.0,
                -0.00041775465696698233,
                -4.362974596862037e-08,
                -1.2151022357848552e-11,
                -4.7588881620421004e-15,
                -2.226614930167366e-18,
                -1.1627237498131586e-21,
            ],
            &mut vec![
                0.0,
                -0.0008355098973052918,
                -1.7444619952659748e-07,
                -7.286557795511902e-11,
                -3.80472772706481e-14,
                -2.2251271876594078e-17,
                1.2789961247944744e-20,
            ],
            &mut vec![
                0.0,
                0.00020861391868413911,
                4.3547247296823945e-08,
                1.515432276542012e-11,
                6.645637323698485e-15,
                3.3399223952510497e-18,
            ],
        );
        assert_eq!(res2.0, 6.046459990680098e-17);
        assert_eq!(res2.1, 0.7094236375834774);
        assert_eq!(res2.2, 0.7047823085448635);
        assert_eq!(res2.3, 0.024682338906797385);
        assert_eq!(res2.4, -0.02467983282954624);
        assert_eq!(res2.5, 0.9996954065371639);
        assert_eq!(res2.6, -0.0);
        assert_eq!(res2.7, 1.0);
        assert_eq!(res2.8, 0.0008355096040059597);
        assert_eq!(res2.9, -5.870849152149326e-05);
        assert_eq!(res2.10, 0.03490027216297455);
    }

    #[test]
    fn test_inverse_start() {
        let geod = Geodesic::new(WGS84_A, WGS84_F);
        let res = geod._InverseStart(
            -0.017393909556108908,
            0.9998487145115275,
            1.0000010195104125,
            -0.0,
            1.0,
            1.0,
            0.017453292519943295,
            0.01745240643728351,
            0.9998476951563913,
            &mut vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            &mut vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        );
        assert_eq!(res.0, -1.0);
        assert_eq!(res.1, 0.7095310092765433);
        assert_eq!(res.2, 0.7046742132893822);
        assert_eq!(res.3.is_nan(), true);
        assert_eq!(res.4.is_nan(), true);
        assert_eq!(res.5, 1.0000002548969817);
    }

    #[test]
    fn test_lengths() {
        // Results taken from the python implementation
        let geod = Geodesic::new(WGS84_A, WGS84_F);
        let mut c1a = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let mut c2a = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let res1 = geod._Lengths(
            0.0008355095326524276,
            0.024682339962725352,
            -0.024679833885152578,
            0.9996954065111039,
            1.0000010195104125,
            -0.0,
            1.0,
            1.0,
            0.9998487145115275,
            1.0,
            4101,
            &mut c1a,
            &mut c2a,
        );
        assert_eq!(res1.0.is_nan(), true);
        assert_eq!(res1.1, 0.024679842274314294);
        assert_eq!(res1.2, 0.0016717180169067588);
        assert_eq!(res1.3.is_nan(), true);
        assert_eq!(res1.4.is_nan(), true);

        let res2 = geod._Lengths(
            0.0008355096040059597,
            0.024682338906797385,
            -0.02467983282954624,
            0.9996954065371639,
            1.0000010195104125,
            -0.0,
            1.0,
            1.0,
            0.9998487145115275,
            1.0,
            4101,
            &mut vec![
                0.0,
                -0.00041775465696698233,
                -4.362974596862037e-08,
                -1.2151022357848552e-11,
                -4.7588881620421004e-15,
                -2.226614930167366e-18,
                -1.1627237498131586e-21,
            ],
            &mut vec![
                0.0,
                -0.0008355098973052918,
                -1.7444619952659748e-07,
                -7.286557795511902e-11,
                -3.80472772706481e-14,
                -2.2251271876594078e-17,
                1.2789961247944744e-20,
            ],
        );
        assert_eq!(res2.0.is_nan(), true);
        assert_eq!(res2.1, 0.02467984121870759);
        assert_eq!(res2.2, 0.0016717181597332804);
        assert_eq!(res2.3.is_nan(), true);
        assert_eq!(res2.4.is_nan(), true);

        let res3 = geod._Lengths(
            0.0008355096040059597,
            0.024682338906797385,
            -0.02467983282954624,
            0.9996954065371639,
            1.0000010195104125,
            -0.0,
            1.0,
            1.0,
            0.9998487145115275,
            1.0,
            1920,
            &mut vec![
                0.0,
                -0.00041775469264372037,
                -4.362975342068502e-08,
                -1.215102547098435e-11,
                -4.758889787701359e-15,
                -2.2266158809456692e-18,
                -1.1627243456014359e-21,
            ],
            &mut vec![
                0.0,
                -0.0008355099686589174,
                -1.744462293162189e-07,
                -7.286559662008413e-11,
                -3.804729026574989e-14,
                -2.2251281376754273e-17,
                1.2789967801615795e-20,
            ],
        );
        assert_eq!(res3.0, 0.024682347295447677);
        assert_eq!(res3.1.is_nan(), true);
        assert_eq!(res3.2.is_nan(), true);
        assert_eq!(res3.3.is_nan(), true);
        assert_eq!(res3.4.is_nan(), true);
    }

    #[test]
    fn test_goed__C4f() {
        let geod = Geodesic::new(WGS84_A, WGS84_F);
        let mut c = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        geod._C4f(0.12, &mut c);
        assert_eq!(
            c,
            vec![
                0.6420952961066771,
                0.0023680700061156517,
                9.96704067834604e-05,
                5.778187189466089e-06,
                3.9979026199316593e-07,
                3.2140078103714466e-08,
                7.0
            ]
        );
    }

    #[test]
    fn test_goed__C3f() {
        let geod = Geodesic::new(WGS84_A, WGS84_F);
        let mut c = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        geod._C3f(0.12, &mut c);

        assert_eq!(
            c,
            vec![
                1.0,
                0.031839442894193756,
                0.0009839921354137713,
                5.0055242248766214e-05,
                3.1656788204092044e-06,
                2.0412e-07,
                7.0
            ]
        );
    }

    #[test]
    fn test_goed__A3f() {
        let geod = Geodesic::new(WGS84_A, WGS84_F);
        assert_eq!(geod._A3f(0.12), 0.9363788874000158);
    }

    #[test]
    fn test_geod_init() {
        // Check that after the init the variables are correctly set.
        // Actual values are taken from the python implementation
        let geod = Geodesic::new(WGS84_A, WGS84_F);
        assert_eq!(geod.a, 6378137.0, "geod.a wrong");
        assert_eq!(geod.f, 0.0033528106647474805, "geod.f wrong");
        assert_eq!(geod._f1, 0.9966471893352525, "geod._f1 wrong");
        assert_eq!(geod._e2, 0.0066943799901413165, "geod._e2 wrong");
        assert_eq!(geod._ep2, 0.006739496742276434, "geod._ep2 wrong");
        assert_eq!(geod._n, 0.0016792203863837047, "geod._n wrong");
        assert_eq!(geod._b, 6356752.314245179, "geod._b wrong");
        assert_eq!(geod._c2, 40589732499314.76, "geod._c2 wrong");
        assert_eq!(geod._etol2, 3.6424611488788524e-08, "geod._etol2 wrong");
        assert_eq!(
            geod._A3x,
            vec![
                -0.0234375,
                -0.046927475637074494,
                -0.06281503005876607,
                -0.2502088451303832,
                -0.49916038980680816,
                1.0
            ],
            "geod._A3x wrong"
        );

        assert_eq!(
            geod._C3x,
            vec![
                0.0234375,
                0.03908873781853724,
                0.04695366939653196,
                0.12499964752736174,
                0.24958019490340408,
                0.01953125,
                0.02345061890926862,
                0.046822392185686165,
                0.062342661206936094,
                0.013671875,
                0.023393770302437927,
                0.025963026642854565,
                0.013671875,
                0.01362595881755982,
                0.008203125
            ],
            "geod._C3x wrong"
        );
        assert_eq!(
            geod._C4x,
            vec![
                0.00646020646020646,
                0.0035037627212872787,
                0.034742279454780166,
                -0.01921732223244865,
                -0.19923321555984239,
                0.6662190894642603,
                0.000111000111000111,
                0.003426620602971002,
                -0.009510765372597735,
                -0.01893413691235592,
                0.0221370239510936,
                0.0007459207459207459,
                -0.004142006291321442,
                -0.00504225176309005,
                0.007584982177746079,
                -0.0021565735851450138,
                -0.001962613370670692,
                0.0036104265913438913,
                -0.0009472009472009472,
                0.0020416649913317735,
                0.0012916376552740189
            ],
            "geod._C4x wrong"
        );
    }
}
