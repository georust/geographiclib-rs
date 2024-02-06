## 0.2.4

* Performance improvements for direct and indirect geodesic calculations
* Remove lazy_static dependency
* Thanks to contributor @valarauca for their work on performance and code cleanup.

## 0.2.3

### New Features

* Added `PolygonArea` to allow calculating perimeter and area of a polygon on a geodesic.
* Added `accurate` feature (enabled by default) for highly accurate `PolygonArea` calculations.

## 0.2.2

This patch release includes only small clippy/lint fixes, and should not
appreciably affect behavior. Thanks to @Quba1 for the cleanup!

- <https://github.com/georust/geographiclib-rs/pull/45>

## 0.2.1

This patch release has no API changes, but fixes some correctness issues around
some of the edge cases represented in Karney's GeodSolveXX tests.

My gratitude to Charles Karney for publishing good tests and especially to
Stony Lohr (@stonylohr) for implementing them here.

And also thank you to Stony for building out our CI and your other project
maintenance work!

<3

### Fixed

* Correctness fixes
  - <https://github.com/georust/geographiclib-rs/pull/34>
  - <https://github.com/georust/geographiclib-rs/pull/35>
  - <https://github.com/georust/geographiclib-rs/pull/37>
  - <https://github.com/georust/geographiclib-rs/pull/41>

### Internal Work

* Set up CI to run against Karney's GeodTest data
  - <https://github.com/georust/geographiclib-rs/pull/17>
  - <https://github.com/georust/geographiclib-rs/pull/18>

* Fix bug in tests
  - <https://github.com/georust/geographiclib-rs/pull/33>
  - <https://github.com/georust/geographiclib-rs/pull/38>

* Code cleanup
  - <https://github.com/georust/geographiclib-rs/pull/39>

## 0.2.0

* Initial release

