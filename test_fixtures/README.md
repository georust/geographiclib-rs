# Test fixtures

## Geodesic

### `GeodTest-100.dat` (built-in)

Small built-in subset of Karney's geodesic test data. The original datasets are published at: <https://sourceforge.net/projects/geographiclib/files/testdata/>

GeographicLib provides the full `GeodTest.dat` dataset and `GeodTest-short.dat`, which samples every 50th row from the full dataset.

`GeodTest-100.dat` is a random 100-row sample from `GeodTest-short.dat`.

### Full and short fixtures

Use `script/download-test-data.sh` and the `test_full` or `test_short` feature flags to run against the original `GeodTest.dat` or `GeodTest-short.dat`.


## Geocentric

### `geocentric_wgs84.dat` (built-in)

GeographicLib does not publish a dedicated Geocentric test dataset, so this file is generated directly from the C++ implementation.

To regenerate:

```sh
./script/generate-geocentric-fixture.sh
```

Columns:

```text
input_lat input_lon input_height x y z output_lat output_lon output_height
```
