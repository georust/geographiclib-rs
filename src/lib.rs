mod geodesic;
pub use geodesic::Geodesic;

pub mod geodesiccapability;
pub use geodesiccapability as capability;

mod geodesicline;
mod geomath;

// TODO only use this for benchmark
use geographiclib;
