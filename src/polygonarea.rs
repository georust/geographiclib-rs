use crate::geodesic::DirectGeodesic;
use crate::geodesic::InverseGeodesic;
use crate::geomath::ang_diff;
use crate::geomath::ang_normalize;
use crate::Geodesic;

#[cfg(feature = "accurate")]
use accurate::traits::*;

/// Clockwise or CounterClockise winding
#[derive(Debug, Clone)]
pub enum Winding {
    Clockwise,
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
////
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
            let (s12, _azi1, _azi2, _m12, _M12, _M21, S12, _a12) =
                self.geoid
                    .inverse(self.latest_lat, self.latest_lon, lat, lon);
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
        let (lat, lon, _azi2, _m12, _M12, _M21, S12, _a12) =
            self.geoid
                .direct(self.latest_lat, self.latest_lon, azimuth, distance);
        self.perimetersum += distance;
        self.areasum += S12;
        self.crossings += PolygonArea::transit(self.latest_lon, lon);
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
    /// let (_perimeter, mut area, _count) = pa.compute(false);
    ///
    /// // Over 5 trillion square meters!
    /// assert_eq!(area, 510053312945726.94); 
    /// ```
    pub fn compute(mut self, sign: bool) -> (f64, f64, usize) {
        #[allow(non_snake_case)]
        let (s12, _azi1, _azi2, _m12, _M12, _M21, S12, _a12) = self.geoid.inverse(
            self.latest_lat,
            self.latest_lon,
            self.initial_lat,
            self.initial_lon,
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

        return (perimetersum, areasum, self.num);
    }

    // Return 1 or -1 if crossing prime meridian in east or west direction.
    // Otherwise return zero.  longitude = +/-0 considered to be positive.
    fn transit(lon1: f64, lon2: f64) -> i64 {
        let (lon12, _lon12s) = ang_diff(lon1, lon2);
        let lon1 = ang_normalize(lon1);
        let lon2 = ang_normalize(lon2);

        // Translation from the following cpp code:
        //  https://github.com/geographiclib/geographiclib/blob/main/src/PolygonArea.cpp#L22
        //
        //  lon12 > 0 && ((lon1 < 0 && lon2 >= 0) ||
        //  (lon1 > 0 && lon2 == 0)) ? 1 :
        //  (lon12 < 0 && lon1 >= 0 && lon2 < 0 ? -1 : 0);

        if lon12 > 0.0 && ((lon1 < 0.0 && lon2 >= 0.0) || (lon1 > 0.0 && lon2 == 0.0)) {
            return 1;
        } else if lon12 < 0.0 && lon1 >= 0.0 && lon2 < 0.0 {
            return -1;
        } else {
            return 0;
        }
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
        }
        else {
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
    fn test_add_edge() {
        let geoid = Geodesic::wgs84();
        let mut pa = PolygonArea::new(&geoid, Winding::CounterClockwise);

        pa.add_point(0.0, 0.0);

        let (s12, azi1, _, _) = geoid.inverse(0.0, 0.0, 0.0, 1.0);
        pa.add_edge(azi1, s12);

        let (s12, azi1, _, _) = geoid.inverse(0.0, 1.0, 1.0, 1.0);
        pa.add_edge(azi1, s12);

        let (s12, azi1, _, _) = geoid.inverse(1.0, 1.0, 1.0, 0.0);
        pa.add_edge(azi1, s12);

        let (perimeter, area, _count) = pa.compute(true);

        assert_relative_eq!(perimeter, 443770.917, epsilon = 1.0e-3);
        assert_relative_eq!(area, 12308778361.469, epsilon = 1.0e-3);
    }

    #[test]
    fn test_planimeter0() {
        // Copied from https://github.com/geographiclib/geographiclib-octave/blob/main/inst/geographiclib_test.m#L644

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
    fn test_planimeter0_add_edge() {
        // Same test as above, but using add_edge

        let geoid = Geodesic::wgs84();

        let mut pa = PolygonArea::new(&geoid, Winding::CounterClockwise);
        pa.add_point(89.0, 0.0);
        let (s12, azi1, _, _) = geoid.inverse(89.0, 0.0, 89.0, 90.0);
        pa.add_edge(azi1, s12);
        let (s12, azi1, _, _) = geoid.inverse(89.0, 90.0, 89.0, 180.0);
        pa.add_edge(azi1, s12);
        let (s12, azi1, _, _) = geoid.inverse(89.0, 180.0, 89.0, 270.0);
        pa.add_edge(azi1, s12);
        let (perimeter, area, _) = pa.compute(true);
        assert_relative_eq!(perimeter, 631819.8745, epsilon = 1.0e-4);
        assert_relative_eq!(area, 24952305678.0, epsilon = 1.0);

        let mut pa = PolygonArea::new(&geoid, Winding::CounterClockwise);
        pa.add_point(-89.0, 0.0);
        let (s12, azi1, _, _) = geoid.inverse(-89.0, 0.0, -89.0, 90.0);
        pa.add_edge(azi1, s12);
        let (s12, azi1, _, _) = geoid.inverse(-89.0, 90.0, -89.0, 180.0);
        pa.add_edge(azi1, s12);
        let (s12, azi1, _, _) = geoid.inverse(-89.0, 180.0, -89.0, 270.0);
        pa.add_edge(azi1, s12);
        let (perimeter, area, _) = pa.compute(true);
        assert_relative_eq!(perimeter, 631819.8745, epsilon = 1.0e-4);
        assert_relative_eq!(area, -24952305678.0, epsilon = 1.0);

        let mut pa = PolygonArea::new(&geoid, Winding::CounterClockwise);
        pa.add_point(0.0, -1.0);
        let (s12, azi1, _, _) = geoid.inverse(0.0, -1.0, -1.0, 0.0);
        pa.add_edge(azi1, s12);
        let (s12, azi1, _, _) = geoid.inverse(-1.0, 0.0, 0.0, 1.0);
        pa.add_edge(azi1, s12);
        let (s12, azi1, _, _) = geoid.inverse(0.0, 1.0, 1.0, 0.0);
        pa.add_edge(azi1, s12);
        let (perimeter, area, _) = pa.compute(true);
        assert_relative_eq!(perimeter, 627598.2731, epsilon = 1.0e-4);
        assert_relative_eq!(area, 24619419146.0, epsilon = 1.0);

        let mut pa = PolygonArea::new(&geoid, Winding::CounterClockwise);
        pa.add_point(90.0, 0.0);
        let (s12, azi1, _, _) = geoid.inverse(90.0, 0.0, 0.0, 0.0);
        pa.add_edge(azi1, s12);
        let (s12, azi1, _, _) = geoid.inverse(0.0, 0.0, 0.0, 90.0);
        pa.add_edge(azi1, s12);
        let (perimeter, area, _) = pa.compute(true);
        assert_relative_eq!(perimeter, 30022685.0, epsilon = 1.0);
        assert_relative_eq!(area, 63758202715511.0, epsilon = 1.0);
    }

    #[test]
    fn test_planimeter5() {
        // Copied from https://github.com/geographiclib/geographiclib-octave/blob/main/inst/geographiclib_test.m#L670

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
        // Copied from https://github.com/geographiclib/geographiclib-octave/blob/main/inst/geographiclib_test.m#L679

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
        // Copied from https://github.com/geographiclib/geographiclib-octave/blob/main/inst/geographiclib_test.m#L701
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
        // Copied from https://github.com/geographiclib/geographiclib-octave/blob/main/inst/geographiclib_test.m#L710
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
        // Copied from https://github.com/geographiclib/geographiclib-octave/blob/main/inst/geographiclib_test.m#L719

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
        // Copied from https://github.com/geographiclib/geographiclib-octave/blob/main/inst/geographiclib_test.m#LL728C14-L728C26

        let geoid = Geodesic::wgs84();

        let mut pa = PolygonArea::new(&geoid, Winding::CounterClockwise);
        pa.add_point(2.0, 1.0);
        pa.add_point(1.0, 2.0);
        pa.add_point(3.0, 3.0);
        let (_, area, _) = pa.compute(true);
        assert_relative_eq!(area, 18454562325.45119);

        // Switching the winding
        let mut pa = PolygonArea::new(&geoid, Winding::Clockwise);
        pa.add_point(2.0, 1.0);
        pa.add_point(1.0, 2.0);
        pa.add_point(3.0, 3.0);
        let (_, area, _) = pa.compute(true);
        assert_relative_eq!(area, -18454562325.45119);

        // Swaping lat and lon
        let mut pa = PolygonArea::new(&geoid, Winding::CounterClockwise);
        pa.add_point(1.0, 2.0);
        pa.add_point(2.0, 1.0);
        pa.add_point(3.0, 3.0);
        let (_, area, _) = pa.compute(true);
        assert_relative_eq!(area, -18454562325.45119);
    }

    #[test]
    fn test_planimeter21() {
        // Copied from https://github.com/geographiclib/geographiclib-octave/blob/main/inst/geographiclib_test.m#L752

        // Testing degenrate polygons.
        let geoid = Geodesic::wgs84();

        let mut pa = PolygonArea::new(&geoid, Winding::CounterClockwise);
        pa.add_point(1.0, 1.0);
        let (perimeter, area, _) = pa.compute(true);
        assert_relative_eq!(perimeter, 0.0);
        assert_relative_eq!(area, 0.0);
    }

    #[test]
    fn test_planimeter29() {
        // Check transitdirect vs transit zero handling consistency

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
