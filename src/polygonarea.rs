use crate::Geodesic;
use crate::geodesic::InverseGeodesic;
use crate::geodesic::DirectGeodesic;
use crate::geomath::ang_normalize;
use crate::geomath::ang_diff;

/// Compute the perimeter and area of a polygon on a Geodesic.
#[derive(Debug, Clone)]
pub struct PolygonArea<'a> {
    geoid: &'a Geodesic,
    num: usize,
    areasum: f64,
    perimetersum: f64,
    crossings: f64,
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
/// use geographiclib_rs::Geodesic;
/// use geographiclib_rs::polygonarea::PolygonArea;
///
/// let g = Geodesic::wgs84();
/// let mut pa = PolygonArea::new(&g);
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
    pub fn new(geoid: &'a Geodesic) -> PolygonArea {
        PolygonArea {
            geoid,
            num: 0,
            areasum: 0.0,
            perimetersum: 0.0,
            crossings: 0.0,
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
            let (s12, _azi1, _azi2, _m12, _M12, _M21, S12, _a12) = self.geoid.inverse(self.latest_lat, self.latest_lon, lat, lon);
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
        let (lat, lon, _azi2, _m12, _M12, _M21, S12, _a12) = self.geoid.direct(self.latest_lat, self.latest_lon, azimuth, distance);
        self.perimetersum += distance;
        self.areasum += S12;
        self.crossings += PolygonArea::transit(self.latest_lon, lon);
        self.latest_lat = lat;
        self.latest_lon = lon;
        self.num += 1;
    }

    /// Return perimeter and area of polygon in meters / meters², consuming PolygonArea.
    pub fn compute(mut self) -> (f64, f64) {
        #[allow(non_snake_case)]
        let (s12, _azi1, _azi2, _m12, _M12, _M21, S12, _a12) = self.geoid.inverse(self.latest_lat, self.latest_lon, self.initial_lat, self.initial_lon);
        self.perimetersum += s12;
        self.areasum += S12;
        self.crossings += PolygonArea::transit(self.latest_lon, self.initial_lon);

        // TODO: Properly take into account crossings when calculating area and perimeter.

        return (self.perimetersum, self.areasum.abs())
    }

    // Return 1 or -1 if crossing prime meridian in east or west direction.
    // Otherwise return zero.  longitude = +/-0 considered to be positive.
    fn transit(lon1: f64, lon2: f64) -> f64 {
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
            return 1.0;
        } else if lon12 < 0.0 && lon1 >= 0.0 && lon2 < 0.0 {
            return -1.0
        }
        else {
            return 0.0;
        }
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
        let mut pa = PolygonArea::new(&geoid);

        pa.add_point(0.0, 0.0);
        pa.add_point(0.0, 1.0);
        pa.add_point(1.0, 1.0);
        pa.add_point(1.0, 0.0);

        let (perimeter, area) = pa.compute();

        assert_relative_eq!(perimeter, 443770.917, epsilon = 1.0e-3);
        assert_relative_eq!(area, 12308778361.469, epsilon = 1.0e-3);
    }
}