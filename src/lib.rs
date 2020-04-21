mod geodesic;
pub use geodesic::Geodesic;

pub mod geodesiccapability;
pub use geodesiccapability as capability;

mod geodesicline;
mod geomath;

#[macro_use]
extern crate lazy_static;
