use crate::geodesic::DirectGeodesic;
use crate::geodesic::InverseGeodesic;
use crate::geomath::ang_diff;
use crate::geomath::ang_normalize;
use crate::Geodesic;

#[cfg(feature = "accurate")]
use accurate::traits::*;

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

// Documentation for PolygonArea
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
/// let (perimeter, area) = pa.compute();
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

    /// Add an adge to the polygon using an azimuth (in degrees) and a distance (in meters). This can only be called after at least one point has been added.
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
        self.crossings += PolygonArea::transit_direct(self.latest_lon, lon);
        self.latest_lat = lat;
        self.latest_lon = lon;
        self.num += 1;
    }

    /// Return perimeter and area of polygon in meters / meters², consuming PolygonArea.
    pub fn compute(mut self) -> (f64, f64) {
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
        let areasum = self.reduce_area(areasum);

        return (perimetersum, areasum);
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

    fn transit_direct(lon1: f64, lon2: f64) -> i64 {
        // Translation of the following cpp code:

        // lon1 = remainder(lon1, real(2 * Math::td));
        // lon2 = remainder(lon2, real(2 * Math::td));
        // return ( (lon2 >= 0 && lon2 < Math::td ? 0 : 1) -
        //         (lon1 >= 0 && lon1 < Math::td ? 0 : 1) );

        let lon1 = lon1 % 720.0;
        let lon2 = lon2 % 720.0;

        let a = {
            if lon2 >= 0.0 && lon2 < 360.0 {
                0
            } else {
                1
            }
        };

        let b = {
            if lon1 >= 0.0 && lon1 < 360.0 {
                0
            } else {
                1
            }
        };

        return a - b;
    }

    fn reduce_area(&self, area: f64) -> f64 {
        let area0 = self.geoid._c2 * 4.0 * std::f64::consts::PI; // Area of earth
        let mut area = area % area0;

        // Translation of the following cpp code:
        // if (crossings & 1) area += (area < 0 ? 1 : -1) * _area0/2;
        if self.crossings % 2 != 0 {
            if area < 0.0 {
                area += area0 / 2.0;
            } else {
                area -= area0 / 2.0;
            }
        }

        // Area is with the clockwise sense. If needed convert to counter-clockwise convention.
        area = match self.winding {
            Winding::Clockwise => area,
            Winding::CounterClockwise => -area,
        };

        // Put area in [0, area0)
        if area >= area0 {
            area -= area0;
        } else if area < 0.0 {
            area += area0;
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
    fn test_polygonarea() {
        let geoid = Geodesic::wgs84();
        let mut pa = PolygonArea::new(&geoid, Winding::CounterClockwise);

        pa.add_point(0.0, 0.0);
        pa.add_point(0.0, 1.0);
        pa.add_point(1.0, 1.0);
        pa.add_point(1.0, 0.0);

        let (perimeter, area) = pa.compute();

        assert_relative_eq!(perimeter, 443770.917, epsilon = 1.0e-3);
        assert_relative_eq!(area, 12308778361.469, epsilon = 1.0e-3);
    }
}
