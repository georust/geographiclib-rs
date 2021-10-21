mod geodesic;
mod gnomonic;
pub use geodesic::{DirectGeodesic, Geodesic, InverseGeodesic};
pub use gnomonic::Gnomonic;

pub mod geodesiccapability;
pub use geodesiccapability as capability;

mod geodesicline;
mod geomath;

#[macro_use]
extern crate lazy_static;
