# geographiclib-rs

A subset of [geographiclib](https://geographiclib.sourceforge.io/) implemented in Rust.

[Documentation](https://docs.rs/geographiclib-rs)

Currently this implements the direct and the inverse geodesic calculations.

If instead you are looking for Rust bindings to [Karney's C++ implementation](https://sourceforge.net/projects/geographiclib/), see [https://crates.io/geographiclib](https://crates.io/crates/geographiclib).

## Examples

```rust
// Determine the point 10000 km NE of JFK - the "direct" geodesic calculation.
use geographiclib_rs::{Geodesic, DirectGeodesic};

let g = Geodesic::wgs84();
let jfk_lat = 40.64;
let jfk_lon = -73.78;
let northeast_azimuth = 45.0;

let (lat, lon, az) = g.direct(jfk_lat, jfk_lon, northeast_azimuth, 10e6);

use approx::assert_relative_eq;
assert_relative_eq!(lat, 32.621100463725796);
assert_relative_eq!(lon, 49.05248709295982);
assert_relative_eq!(az,  140.4059858768007);
```

```rust
// Determine the distance between two points - the "inverse" geodesic calculation.
use geographiclib_rs::{Geodesic, InverseGeodesic};

let g = Geodesic::wgs84();
let p1 = (34.095925, -118.2884237);
let p2 = (59.4323439, 24.7341649);
let s12: f64 = g.inverse(p1.0, p1.1, p2.0, p2.1);

use approx::assert_relative_eq;
assert_relative_eq!(s12, 9094718.72751138);
``` 

```rust
// Determine the perimeter and area of a polygon.
use geographiclib_rs::{Geodesic, PolygonArea, Winding};

let g = Geodesic::wgs84();
let mut pa = PolygonArea::new(&g, Winding::CounterClockwise);
pa.add_point(0.0, 0.0);
pa.add_point(0.0, 1.0);
pa.add_point(1.0, 1.0);
pa.add_point(1.0, 0.0);

let (perimeter_m, area_m_squared, num_points) = pa.compute(false);

use approx::assert_relative_eq;
assert_relative_eq!(perimeter_m, 443770.91724830196);
assert_relative_eq!(area_m_squared, 12308778361.469452);
assert_eq!(num_points, 4);
``` 

```rust
// Determine the distance between rovers Pathfinder and Curiosity on Mars
use geographiclib_rs::{Geodesic, InverseGeodesic};

let mars = Geodesic::new(3396190.0, 1.0 / 169.8944472);
let pathfinder = (19.26, 326.75);
let curiosity = (-4.765700445, 137.39820983);
let distance_m: f64 = mars.inverse(curiosity.0, curiosity.1, pathfinder.0, pathfinder.1);

assert_eq!(distance_m.round(), 9639113.0);
```

## Features

1. `accurate`: Enabled by default. Use the [`accurate`](https://docs.rs/accurate/latest/accurate/) crate to provide high accuracy polygon areas and perimeters in `PolygonArea`. Can be disabled for better performance or when `PolygonArea` is not being used.

## Benchmarking

To compare the direct and inverse geodesic calculation against the [geographiclib c bindings](https://github.com/savage13/geographiclib), run:

```shell
cargo bench
```

Which produces output like:

```text
direct (c wrapper)/default
                        time:   [24.046 µs 24.071 µs 24.099 µs]

direct (rust impl)/default
                        time:   [26.129 µs 26.168 µs 26.211 µs]

inverse (c wrapper)/default
                        time:   [45.061 µs 45.141 µs 45.227 µs]

inverse (rust impl)/default
                        time:   [67.739 µs 67.796 µs 67.865 µs]
```

Showing that, at least in this benchmark, the Rust implementation is 10-50% slower than the c bindings.
