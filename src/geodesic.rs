use crate::geodesiccapability;
use crate::geomath;

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

    pub fn _C3f(&self, eps: f64, c: Vec<f64>) -> Vec<f64> {
        let mut c = c;
        let mut mult = 1.0;
        let mut o = 0.0;
        for l in 1..self.GEODESIC_ORDER {
            let m = self.GEODESIC_ORDER - l - 1;
            mult *= eps;
            c[l as usize] = mult * geomath::polyval(m, &self._C3x, o as usize, eps);
            o += m as f64 + 1.0;
        }
        c
    }

    pub fn _C4f(&self, eps: f64, c: Vec<f64>) -> Vec<f64> {
        let mut c = c;
        let mut mult = 1.0;
        let mut o = 0.0;
        for l in 0..self.GEODESIC_ORDER {
            let m = self.GEODESIC_ORDER - l - 1;
            c[l as usize] = mult * geomath::polyval(m, &self._C4x, o as usize, eps);
            o += m as f64 + 1.0;
            mult *= eps;
        }
        c
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_goed__C4f() {
        const WGS84_A: f64 = 6378137.0;
        const WGS84_F: f64 = 1.0 / 298.257223563;
        let geod = Geodesic::new(WGS84_A, WGS84_F);
        assert_eq!(
            geod._C4f(0.12, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]),
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
        const WGS84_A: f64 = 6378137.0;
        const WGS84_F: f64 = 1.0 / 298.257223563;
        let geod = Geodesic::new(WGS84_A, WGS84_F);
        assert_eq!(
            geod._C3f(0.12, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]),
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
        const WGS84_A: f64 = 6378137.0;
        const WGS84_F: f64 = 1.0 / 298.257223563;
        let geod = Geodesic::new(WGS84_A, WGS84_F);
        assert_eq!(geod._A3f(0.12), 0.9363788874000158);
    }

    #[test]
    fn test_geod_init() {
        // Check that after the init the variables are correctly set.
        // Actual values are taken from the python implementation
        const WGS84_A: f64 = 6378137.0;
        const WGS84_F: f64 = 1.0 / 298.257223563;
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
