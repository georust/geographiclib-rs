use crate::geomath::ang_diff;
use crate::geomath::ang_normalize;
use crate::Geodesic;

use crate::geodesic_capability as caps;

const POLYGONAREA_MASK: u64 =
    caps::LATITUDE | caps::LONGITUDE | caps::DISTANCE | caps::AREA | caps::LONG_UNROLL;

#[cfg(feature = "accurate")]
use accurate::traits::*;

/// Clockwise or CounterClockwise winding
///
/// The standard winding of a Simple Feature polygon is counter-clockwise. However, if the polygon is a hole, then the winding is clockwise.
/// ESRI Shapefile polygons are opposite, with the outer-ring being clockwise and holes being counter-clockwise.
#[derive(Debug, Default, Copy, Clone)]
pub enum Winding {
    Clockwise,
    #[default]
    CounterClockwise,
}

/// Compute the perimeter and area of a polygon on a Geodesic.
#[derive(Debug, Clone)]
pub struct PolygonArea<'a> {
    geoid: &'a Geodesic,
    winding: Winding,
    num: usize,

    #[cfg(not(feature = "accurate"))]
    areasum: f64,

    #[cfg(feature = "accurate")]
    areasum: accurate::sum::Sum2<f64>,

    #[cfg(not(feature = "accurate"))]
    perimetersum: f64,

    #[cfg(feature = "accurate")]
    perimetersum: accurate::sum::Sum2<f64>,

    crossings: i64,
    initial_lat: f64,
    initial_lon: f64,
    latest_lat: f64,
    latest_lon: f64,
}

/// PolygonArea can be used to compute the perimeter and area of a polygon on a Geodesic.
///
/// # Example
/// ```rust
/// use geographiclib_rs::{Geodesic, PolygonArea, Winding};
///
/// let g = Geodesic::wgs84();
/// let mut pa = PolygonArea::new(&g, Winding::CounterClockwise);
///
/// pa.add_point(0.0, 0.0);
/// pa.add_point(0.0, 1.0);
/// pa.add_point(1.0, 1.0);
/// pa.add_point(1.0, 0.0);
///
/// let (perimeter, area, _num) = pa.compute(false);
///
/// use approx::assert_relative_eq;
/// assert_relative_eq!(perimeter, 443770.917248302);
/// assert_relative_eq!(area, 12308778361.469452);
/// ```
impl<'a> PolygonArea<'a> {
    /// Create a new PolygonArea using a Geodesic.
    pub fn new(geoid: &'a Geodesic, winding: Winding) -> PolygonArea {
        PolygonArea {
            geoid,
            winding,

            num: 0,

            #[cfg(not(feature = "accurate"))]
            areasum: 0.0,

            #[cfg(feature = "accurate")]
            areasum: accurate::sum::Sum2::zero(),

            #[cfg(not(feature = "accurate"))]
            perimetersum: 0.0,

            #[cfg(feature = "accurate")]
            perimetersum: accurate::sum::Sum2::zero(),

            crossings: 0,
            initial_lat: 0.0,
            initial_lon: 0.0,
            latest_lat: 0.0,
            latest_lon: 0.0,
        }
    }

    /// Add a point to the polygon
    pub fn add_point(&mut self, lat: f64, lon: f64) {
        if self.num == 0 {
            self.initial_lat = lat;
            self.initial_lon = lon;
        } else {
            #[allow(non_snake_case)]
            let (_a12, s12, _salp1, _calp1, _salp2, _calp2, _m12, _M12, _M21, S12) = self
                .geoid
                ._gen_inverse(self.latest_lat, self.latest_lon, lat, lon, POLYGONAREA_MASK);
            self.perimetersum += s12;
            self.areasum += S12;
            self.crossings += PolygonArea::transit(self.latest_lon, lon);
        }
        self.latest_lat = lat;
        self.latest_lon = lon;
        self.num += 1;
    }

    /// Add an edge to the polygon using an azimuth (in degrees) and a distance (in meters). This can only be called after at least one point has been added.
    ///
    /// # Panics
    /// Panics if no points have been added yet.
    pub fn add_edge(&mut self, azimuth: f64, distance: f64) {
        if self.num == 0 {
            panic!("PolygonArea::add_edge: No points added yet");
        }

        #[allow(non_snake_case)]
        let (_a12, lat, lon, _azi2, _s12, _m12, _M12, _M21, S12) = self.geoid._gen_direct(
            self.latest_lat,
            self.latest_lon,
            azimuth,
            false,
            distance,
            POLYGONAREA_MASK,
        );
        self.perimetersum += distance;
        self.areasum += S12;
        self.crossings += PolygonArea::transitdirect(self.latest_lon, lon);
        self.latest_lat = lat;
        self.latest_lon = lon;
        self.num += 1;
    }

    /// Consumes the PolygonArea and returns the following tuple:
    ///  - 0: Perimeter in (meters) of the polygon.
    ///  - 1: Area (meters²) of the polygon.
    ///  - 2: Number of points added to the polygon.
    ///
    /// # Parameters
    ///
    /// - `sign`: Whether to allow negative values for the area.
    ///   - `true`: Interpret an inversely wound polygon to be a "negative" area. This will produce incorrect results if the polygon covers over half the geodesic. See "Interpreting negative area values" below.
    ///   - `false`: Always return a positive area. Inversely wound polygons are treated as if they are always wound the same way as the winding specified during creation of the PolygonArea. This is useful if you are dealing with very large polygons that might cover over half the geodesic, or if you are certain of your winding.
    ///
    /// # Interpreting negative area values
    ///
    /// A negative value can mean one of two things:
    /// 1. The winding of the polygon is opposite the winding specified during creation of the PolygonArea.
    /// 2. The polygon is larger than half the planet. In this case, to get the final area of the polygon, add the area of the planet to the result. If you expect to be dealing with polygons of this size pass `signed = false` to `compute()` to get the correct result.
    ///
    /// # Large polygon example
    /// ```rust
    /// use geographiclib_rs::{Geodesic, PolygonArea, Winding};
    /// let g = Geodesic::wgs84();
    ///
    /// // Describe a polygon that covers all of the earth EXCEPT this small square.
    /// // The outside of the polygon is in this square, the inside of the polygon is the rest of the earth.
    /// let mut pa = PolygonArea::new(&g, Winding::CounterClockwise);
    /// pa.add_point(0.0, 0.0);
    /// pa.add_point(1.0, 0.0);
    /// pa.add_point(1.0, 1.0);
    /// pa.add_point(0.0, 1.0);
    ///
    /// let (_perimeter, area, _count) = pa.compute(false);
    ///
    /// // Over 5 trillion square meters!
    /// assert_eq!(area, 510053312945726.94);
    /// ```
    pub fn compute(mut self, sign: bool) -> (f64, f64, usize) {
        #[allow(non_snake_case)]
        let (_a12, s12, _salp1, _calp1, _salp2, _calp2, _m12, _M12, _M21, S12) =
            self.geoid._gen_inverse(
                self.latest_lat,
                self.latest_lon,
                self.initial_lat,
                self.initial_lon,
                POLYGONAREA_MASK,
            );
        self.perimetersum += s12;
        self.areasum += S12;

        let perimetersum;
        let areasum;

        #[cfg(not(feature = "accurate"))]
        {
            perimetersum = self.perimetersum;
            areasum = self.areasum;
        }
        #[cfg(feature = "accurate")]
        {
            perimetersum = self.perimetersum.sum();
            areasum = self.areasum.sum();
        }

        self.crossings += PolygonArea::transit(self.latest_lon, self.initial_lon);

        // Properly take into account crossings when calculating area.
        let areasum = self.reduce_area(areasum, sign);

        (perimetersum, areasum, self.num)
    }

    /// Check what the perimeter and area would be if this point was added to the polygon without actually adding it
    pub fn test_point(&self, lat: f64, lon: f64, sign: bool) -> (f64, f64, usize) {
        let mut pa = self.clone();
        pa.add_point(lat, lon);
        pa.compute(sign)
    }

    /// Check what the perimeter and area would be if this edge was added to the polygon without actually adding it
    pub fn test_edge(&self, azimuth: f64, distance: f64, sign: bool) -> (f64, f64, usize) {
        let mut pa = self.clone();
        pa.add_edge(azimuth, distance);
        pa.compute(sign)
    }

    // Return 1 or -1 if crossing prime meridian in east or west direction.
    // Otherwise return zero.  longitude = +/-0 considered to be positive.
    fn transit(lon1: f64, lon2: f64) -> i64 {
        let (lon12, _lon12s) = ang_diff(lon1, lon2);
        let lon1 = ang_normalize(lon1);
        let lon2 = ang_normalize(lon2);

        // Translation from the following cpp code:
        //  https://github.com/geographiclib/geographiclib/blob/8bc13eb53acdd8bc4fbe4de212d42dbb29779476/src/PolygonArea.cpp#L22
        //
        //  lon12 > 0 && ((lon1 < 0 && lon2 >= 0) ||
        //  (lon1 > 0 && lon2 == 0)) ? 1 :
        //  (lon12 < 0 && lon1 >= 0 && lon2 < 0 ? -1 : 0);

        if lon12 > 0.0 && ((lon1 < 0.0 && lon2 >= 0.0) || (lon1 > 0.0 && lon2 == 0.0)) {
            1
        } else if lon12 < 0.0 && lon1 >= 0.0 && lon2 < 0.0 {
            -1
        } else {
            0
        }
    }

    fn transitdirect(lon1: f64, lon2: f64) -> i64 {
        // We want to compute exactly: floor(lon2 / 360) - floor(lon1 / 360)
        let lon1 = lon1 % 720.0;
        let lon2 = lon2 % 720.0;

        let a = if 0.0 <= lon2 && lon2 < 360.0 { 0 } else { 1 };

        let b = if 0.0 <= lon1 && lon1 < 360.0 { 0 } else { 1 };

        a - b
    }

    fn reduce_area(&self, area: f64, signed: bool) -> f64 {
        let geoid_area = self.geoid.area(); // Area of the planet
        let mut area = area % geoid_area;

        // Translation of the following cpp code:
        // if (crossings & 1) area += (area < 0 ? 1 : -1) * _area0/2;
        if self.crossings % 2 != 0 {
            if area < 0.0 {
                area += geoid_area / 2.0;
            } else {
                area -= geoid_area / 2.0;
            }
        }

        // Area is with the clockwise sense. If needed convert to counter-clockwise convention.
        area = match self.winding {
            Winding::Clockwise => area,
            Winding::CounterClockwise => -area,
        };

        if signed {
            // Put area in (-geoid_area/2, geoid_area/2]
            if area > geoid_area / 2.0 {
                area -= geoid_area;
            } else if area <= -geoid_area / 2.0 {
                area += geoid_area;
            }
        } else {
            // Negative values mean the polygon is larger than half the planet. Correct for this.
            if area < 0.0 {
                area += geoid_area;
            }
        }

        area
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Geodesic;
    use approx::assert_relative_eq;

    #[test]
    fn test_simple_polygonarea() {
        let geoid = Geodesic::wgs84();
        let mut pa = PolygonArea::new(&geoid, Winding::CounterClockwise);

        pa.add_point(0.0, 0.0);
        pa.add_point(0.0, 1.0);
        pa.add_point(1.0, 1.0);
        pa.add_point(1.0, 0.0);

        let (perimeter, area, _count) = pa.compute(true);

        assert_relative_eq!(perimeter, 443770.917, epsilon = 1.0e-3);
        assert_relative_eq!(area, 12308778361.469, epsilon = 1.0e-3);

        let mut pa = PolygonArea::new(&geoid, Winding::Clockwise);

        pa.add_point(0.0, 0.0);
        pa.add_point(0.0, 1.0);
        pa.add_point(1.0, 1.0);
        pa.add_point(1.0, 0.0);

        let (perimeter, area, _count) = pa.compute(true);

        assert_relative_eq!(perimeter, 443770.917, epsilon = 1.0e-3);
        assert_relative_eq!(area, -12308778361.469, epsilon = 1.0e-3);
    }

    #[test]
    fn test_planimeter0() {
        // Copied from https://github.com/geographiclib/geographiclib-octave/blob/0662e05a432a040a60ab27c779fa09b554177ba9/inst/geographiclib_test.m#L644

        let geoid = Geodesic::wgs84();

        let mut pa = PolygonArea::new(&geoid, Winding::CounterClockwise);
        pa.add_point(89.0, 0.0);
        pa.add_point(89.0, 90.0);
        pa.add_point(89.0, 180.0);
        pa.add_point(89.0, 270.0);
        let (perimeter, area, _count) = pa.compute(true);
        assert_relative_eq!(perimeter, 631819.8745, epsilon = 1.0e-4);
        assert_relative_eq!(area, 24952305678.0, epsilon = 1.0);

        let mut pa = PolygonArea::new(&geoid, Winding::CounterClockwise);
        pa.add_point(-89.0, 0.0);
        pa.add_point(-89.0, 90.0);
        pa.add_point(-89.0, 180.0);
        pa.add_point(-89.0, 270.0);
        let (perimeter, area, _count) = pa.compute(true);
        assert_relative_eq!(perimeter, 631819.8745, epsilon = 1.0e-4);
        assert_relative_eq!(area, -24952305678.0, epsilon = 1.0);

        let mut pa = PolygonArea::new(&geoid, Winding::CounterClockwise);
        pa.add_point(0.0, -1.0);
        pa.add_point(-1.0, 0.0);
        pa.add_point(0.0, 1.0);
        pa.add_point(1.0, 0.0);
        let (perimeter, area, _count) = pa.compute(true);
        assert_relative_eq!(perimeter, 627598.2731, epsilon = 1.0e-4);
        assert_relative_eq!(area, 24619419146.0, epsilon = 1.0);

        let mut pa = PolygonArea::new(&geoid, Winding::CounterClockwise);
        pa.add_point(90.0, 0.0);
        pa.add_point(0.0, 0.0);
        pa.add_point(0.0, 90.0);
        let (perimeter, area, _count) = pa.compute(true);
        assert_relative_eq!(perimeter, 30022685.0, epsilon = 1.0);
        assert_relative_eq!(area, 63758202715511.0, epsilon = 1.0);
    }

    #[test]
    fn test_planimeter5() {
        // Copied from https://github.com/geographiclib/geographiclib-octave/blob/0662e05a432a040a60ab27c779fa09b554177ba9/inst/geographiclib_test.m#L670

        let geoid = Geodesic::wgs84();
        let mut pa = PolygonArea::new(&geoid, Winding::CounterClockwise);
        pa.add_point(89.0, 0.1);
        pa.add_point(89.0, 90.1);
        pa.add_point(89.0, -179.9);
        let (perimeter, area, _) = pa.compute(true);
        assert_relative_eq!(perimeter, 539297.0, epsilon = 1.0);
        assert_relative_eq!(area, 12476152838.5, epsilon = 1.0);
    }

    #[test]
    fn test_planimeter6() {
        // Copied from https://github.com/geographiclib/geographiclib-octave/blob/0662e05a432a040a60ab27c779fa09b554177ba9/inst/geographiclib_test.m#L679

        let geoid = Geodesic::wgs84();

        let mut pa = PolygonArea::new(&geoid, Winding::CounterClockwise);
        pa.add_point(9.0, -0.00000000000001);
        pa.add_point(9.0, 180.0);
        pa.add_point(9.0, 0.0);
        let (perimeter, area, _) = pa.compute(true);
        assert_relative_eq!(perimeter, 36026861.0, epsilon = 1.0);
        assert_relative_eq!(area, 0.0, epsilon = 1.0);

        let mut pa = PolygonArea::new(&geoid, Winding::CounterClockwise);
        pa.add_point(9.0, 0.00000000000001);
        pa.add_point(9.0, 0.0);
        pa.add_point(9.0, 180.0);
        let (perimeter, area, _) = pa.compute(true);
        assert_relative_eq!(perimeter, 36026861.0, epsilon = 1.0);
        assert_relative_eq!(area, 0.0, epsilon = 1.0);

        let mut pa = PolygonArea::new(&geoid, Winding::CounterClockwise);
        pa.add_point(9.0, 0.00000000000001);
        pa.add_point(9.0, 180.0);
        pa.add_point(9.0, 0.0);
        let (perimeter, area, _) = pa.compute(true);
        assert_relative_eq!(perimeter, 36026861.0, epsilon = 1.0);
        assert_relative_eq!(area, 0.0, epsilon = 1.0);

        let mut pa = PolygonArea::new(&geoid, Winding::CounterClockwise);
        pa.add_point(9.0, -0.00000000000001);
        pa.add_point(9.0, 0.0);
        pa.add_point(9.0, 180.0);
        let (perimeter, area, _) = pa.compute(true);
        assert_relative_eq!(perimeter, 36026861.0, epsilon = 1.0);
        assert_relative_eq!(area, 0.0, epsilon = 1.0);
    }

    #[test]
    fn test_planimeter12() {
        // Copied from https://github.com/geographiclib/geographiclib-octave/blob/0662e05a432a040a60ab27c779fa09b554177ba9/inst/geographiclib_test.m#L701
        let geoid = Geodesic::wgs84();

        let mut pa = PolygonArea::new(&geoid, Winding::CounterClockwise);
        pa.add_point(66.562222222, 0.0);
        pa.add_point(66.562222222, 180.0);
        pa.add_point(66.562222222, 360.0);
        let (perimeter, area, _) = pa.compute(true);
        assert_relative_eq!(perimeter, 10465729.0, epsilon = 1.0);
        assert_relative_eq!(area, 0.0);
    }

    #[test]
    fn test_planimeter12r() {
        // Copied from https://github.com/geographiclib/geographiclib-octave/blob/0662e05a432a040a60ab27c779fa09b554177ba9/inst/geographiclib_test.m#L710
        let geoid = Geodesic::wgs84();

        let mut pa = PolygonArea::new(&geoid, Winding::CounterClockwise);

        pa.add_point(66.562222222, -0.0);
        pa.add_point(66.562222222, -180.0);
        pa.add_point(66.562222222, -360.0);

        let (perimeter, area, _) = pa.compute(true);
        assert_relative_eq!(perimeter, 10465729.0, epsilon = 1.0);
        assert_relative_eq!(area, 0.0);
    }

    #[test]
    fn test_planimeter13() {
        // Copied from https://github.com/geographiclib/geographiclib-octave/blob/0662e05a432a040a60ab27c779fa09b554177ba9/inst/geographiclib_test.m#L719

        let geoid = Geodesic::wgs84();

        let mut pa = PolygonArea::new(&geoid, Winding::CounterClockwise);
        pa.add_point(89.0, -360.0);
        pa.add_point(89.0, -240.0);
        pa.add_point(89.0, -120.0);
        pa.add_point(89.0, 0.0);
        pa.add_point(89.0, 120.0);
        pa.add_point(89.0, 240.0);
        let (perimeter, area, _) = pa.compute(true);
        assert_relative_eq!(perimeter, 1160741.0, epsilon = 1.0);
        assert_relative_eq!(area, 32415230256.0, epsilon = 1.0);
    }

    #[test]
    fn test_planimeter15() {
        // Copied from https://github.com/geographiclib/geographiclib-octave/blob/0662e05a432a040a60ab27c779fa09b554177ba9/inst/geographiclib_test.m#LL728C14-L728C26

        let geoid = Geodesic::wgs84();

        let mut pa = PolygonArea::new(&geoid, Winding::CounterClockwise);
        pa.add_point(2.0, 1.0);
        pa.add_point(1.0, 2.0);
        pa.add_point(3.0, 3.0);
        let (_, area, _) = pa.compute(true);
        assert_relative_eq!(area, 18454562325.45119, epsilon = 1.0e-4);

        // Switching the winding
        let mut pa = PolygonArea::new(&geoid, Winding::Clockwise);
        pa.add_point(2.0, 1.0);
        pa.add_point(1.0, 2.0);
        pa.add_point(3.0, 3.0);
        let (_, area, _) = pa.compute(true);
        assert_relative_eq!(area, -18454562325.45119, epsilon = 1.0e-4);

        // Swaping lat and lon
        let mut pa = PolygonArea::new(&geoid, Winding::CounterClockwise);
        pa.add_point(1.0, 2.0);
        pa.add_point(2.0, 1.0);
        pa.add_point(3.0, 3.0);
        let (_, area, _) = pa.compute(true);
        assert_relative_eq!(area, -18454562325.45119, epsilon = 1.0e-4);
    }

    #[test]
    fn test_planimeter19() {
        // Test degenerate polygons

        // Copied from https://github.com/geographiclib/geographiclib-js/blob/57137fdcf4ba56718c64b909b00331754b6efceb/geodesic/test/geodesictest.js#L801
        // Testing degenrate polygons.
        let geoid = Geodesic::wgs84();

        let mut pa = PolygonArea::new(&geoid, Winding::CounterClockwise);
        let (perimeter, area, _) = pa.clone().compute(true);
        assert_relative_eq!(perimeter, 0.0);
        assert_relative_eq!(area, 0.0);

        let (perimeter, area, _) = pa.test_point(1.0, 1.0, true);
        assert_relative_eq!(perimeter, 0.0);
        assert_relative_eq!(area, 0.0);

        let result = std::panic::catch_unwind(|| {
            let (_, _, _) = pa.test_edge(90.0, 1000.0, true);
        });
        assert!(result.is_err());

        pa.add_point(1.0, 1.0);
        let (perimeter, area, _) = pa.clone().compute(true);
        assert_relative_eq!(perimeter, 0.0);
        assert_relative_eq!(area, 0.0);
    }

    #[test]
    fn test_planimeter21() {
        // Copied From: https://github.com/geographiclib/geographiclib-js/blob/57137fdcf4ba56718c64b909b00331754b6efceb/geodesic/test/geodesictest.js#LL836C13-L836C13

        let g = Geodesic::wgs84();

        let lat = 45.0;
        let azi = 39.2144607176828184218;
        let s = 8420705.40957178156285;
        let r = 39433884866571.4277; // Area for one circuit
        let a0 = g.area(); // Ellipsoid area

        let mut pa_clockwise = PolygonArea::new(&g, Winding::Clockwise);
        let mut pa_counter = PolygonArea::new(&g, Winding::CounterClockwise);

        pa_clockwise.add_point(lat, 60.0);
        pa_clockwise.add_point(lat, 180.0);
        pa_clockwise.add_point(lat, -60.0);
        pa_clockwise.add_point(lat, 60.0);
        pa_clockwise.add_point(lat, 180.0);
        pa_clockwise.add_point(lat, -60.0);

        pa_counter.add_point(lat, 60.0);
        pa_counter.add_point(lat, 180.0);
        pa_counter.add_point(lat, -60.0);
        pa_counter.add_point(lat, 60.0);
        pa_counter.add_point(lat, 180.0);
        pa_counter.add_point(lat, -60.0);

        for i in 3..=4 {
            pa_clockwise.add_point(lat, 60.0);
            pa_clockwise.add_point(lat, 180.0);

            pa_counter.add_point(lat, 60.0);
            pa_counter.add_point(lat, 180.0);

            let (_, area, _) = pa_counter.test_point(lat, -60.0, true);
            assert_relative_eq!(area, i as f64 * r, epsilon = 0.5);

            let (_, area, _) = pa_counter.test_point(lat, -60.0, false);
            assert_relative_eq!(area, i as f64 * r, epsilon = 0.5);

            let (_, area, _) = pa_clockwise.test_point(lat, -60.0, true);
            assert_relative_eq!(area, -i as f64 * r, epsilon = 0.5);

            let (_, area, _) = pa_clockwise.test_point(lat, -60.0, false);
            assert_relative_eq!(area, (-i as f64 * r) + a0, epsilon = 0.5);

            let (_, area, _) = pa_counter.test_edge(azi, s, true);
            assert_relative_eq!(area, i as f64 * r, epsilon = 0.5);

            let (_, area, _) = pa_counter.test_edge(azi, s, false);
            assert_relative_eq!(area, i as f64 * r, epsilon = 0.5);

            let (_, area, _) = pa_clockwise.test_edge(azi, s, true);
            assert_relative_eq!(area, -i as f64 * r, epsilon = 0.5);

            let (_, area, _) = pa_clockwise.test_edge(azi, s, false);
            assert_relative_eq!(area, (-i as f64 * r) + a0, epsilon = 0.5);

            pa_clockwise.add_point(lat, -60.0);
            pa_counter.add_point(lat, -60.0);

            let (_, area, _) = pa_counter.clone().compute(true);
            assert_relative_eq!(area, r * i as f64, epsilon = 0.5);

            let (_, area, _) = pa_counter.clone().compute(false);
            assert_relative_eq!(area, r * i as f64, epsilon = 0.5);

            let (_, area, _) = pa_clockwise.clone().compute(true);
            assert_relative_eq!(area, -(r * i as f64), epsilon = 0.5);

            let (_, area, _) = pa_clockwise.clone().compute(false);
            assert_relative_eq!(area, -(r * i as f64) + a0, epsilon = 0.5);
        }
    }

    #[test]
    fn test_planimeter29() {
        // Check transitdirect vs transit zero handling consistency
        // Copied from: https://github.com/geographiclib/geographiclib-js/blob/57137fdcf4ba56718c64b909b00331754b6efceb/geodesic/test/geodesictest.js#L883

        let geoid = Geodesic::wgs84();
        let mut pa = PolygonArea::new(&geoid, Winding::CounterClockwise);
        pa.add_point(0.0, 0.0);
        pa.add_edge(90.0, 1000.0);
        pa.add_edge(0.0, 1000.0);
        pa.add_edge(-90.0, 1000.0);
        let (_, area, _) = pa.compute(true);
        assert_relative_eq!(area, 1000000.0, epsilon = 0.01);
    }
}
