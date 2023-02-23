mod geodesic;
pub use geodesic::{DirectGeodesic, Geodesic, InverseGeodesic};

pub mod geodesiccapability;
pub use geodesiccapability as capability;

mod geodesicline;
mod geomath;
pub mod polygonarea;

#[macro_use]
extern crate lazy_static;
