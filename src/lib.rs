//! A subset of [geographiclib](https://geographiclib.sourceforge.io/) implemented in Rust.
//! 
//! This library is a subset of the [geographiclib](https://geographiclib.sourceforge.io/) library
//! 
//! # Examples
//! 
//! ```rust
//! // Determine the point 10000 km NE of JFK - the "direct" geodesic calculation.
//! use geographiclib_rs::{Geodesic, DirectGeodesic};
//! 
//! let g = Geodesic::wgs84();
//! let jfk_lat = 40.64;
//! let jfk_lon = -73.78;
//! let northeast_azimuth = 45.0;
//! 
//! let (lat, lon, az) = g.direct(jfk_lat, jfk_lon, northeast_azimuth, 10e6);
//! 
//! use approx::assert_relative_eq;
//! assert_relative_eq!(lat, 32.621100463725796);
//! assert_relative_eq!(lon, 49.05248709295982);
//! assert_relative_eq!(az,  140.4059858768007);
//! ```
//! 
//! ```rust
//! // Determine the distance between two points - the "inverse" geodesic calculation.
//! use geographiclib_rs::{Geodesic, InverseGeodesic};
//! 
//! let g = Geodesic::wgs84();
//! let p1 = (34.095925, -118.2884237);
//! let p2 = (59.4323439, 24.7341649);
//! let s12: f64 = g.inverse(p1.0, p1.1, p2.0, p2.1);
//! 
//! use approx::assert_relative_eq;
//! assert_relative_eq!(s12, 9094718.72751138);
//! ``` 

mod geodesic;
pub use geodesic::{DirectGeodesic, Geodesic, InverseGeodesic};

pub mod geodesiccapability;
pub use geodesiccapability as capability;

mod geodesicline;
mod geomath;
mod polygonarea;
pub use polygonarea::PolygonArea;
pub use polygonarea::Winding;

#[macro_use]
extern crate lazy_static;
