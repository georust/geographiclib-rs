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

let (perimeter_m, area_m_squared) = pa.compute();

use approx::assert_relative_eq;
assert_relative_eq!(perimeter_m, 443770.91724830196);
assert_relative_eq!(area_m_squared, 12308778361.469452);
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
     Running target/release/deps/geodesic_benchmark-af6ba4f7be913514
direct (c wrapper)      time:   [34.852 us 34.937 us 35.023 us]                                
                        change: [+1.1137% +1.7246% +2.2864%] (p = 0.00 < 0.05)
                        Performance has regressed.
Found 9 outliers among 100 measurements (9.00%)
  3 (3.00%) low mild
  6 (6.00%) high mild

direct (rust impl)      time:   [48.862 us 48.959 us 49.059 us]                                
                        change: [+0.0149% +0.8003% +1.5464%] (p = 0.04 < 0.05)
                        Change within noise threshold.
Found 9 outliers among 100 measurements (9.00%)
  1 (1.00%) low mild
  4 (4.00%) high mild
  4 (4.00%) high severe

inverse (c wrapper)     time:   [70.875 us 71.138 us 71.464 us]                                
                        change: [+0.6259% +1.1321% +1.6653%] (p = 0.00 < 0.05)
                        Change within noise threshold.
Found 8 outliers among 100 measurements (8.00%)
  1 (1.00%) high mild
  7 (7.00%) high severe

inverse (rust impl)     time:   [103.66 us 104.07 us 104.58 us]                                
                        change: [-1.0415% -0.0086% +1.0291%] (p = 0.99 > 0.05)
                        No change in performance detected.
Found 7 outliers among 100 measurements (7.00%)
  1 (1.00%) low mild
  6 (6.00%) high severe
```

Showing that, at least in this benchmark, the Rust implementation is 40-50% slower than the c bindings.
