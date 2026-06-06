# Test fixtures

## Geodesic

### `GeodTest-100.dat` (built-in)

Small built-in subset of Karney's geodesic test data. The original datasets are published at: <https://sourceforge.net/projects/geographiclib/files/testdata/>

GeographicLib provides the full `GeodTest.dat` dataset and `GeodTest-short.dat`, which samples every 50th row from the full dataset.

`GeodTest-100.dat` is a random 100-row sample from `GeodTest-short.dat`.

### Full and short fixtures

Use `script/download-test-data.sh` and the `test_full` or `test_short` feature flags to run against the original `GeodTest.dat` or `GeodTest-short.dat`.


## Transverse Mercator

### `TMcoords-excerpt.dat` (built-in)

Subset of Karney's [`TMcoords.dat`](https://geographiclib.sourceforge.io/C++/doc/transversemercator.html#testmerc), sampled by `script/excerpt-tmcoords.py`.  Rows that require `TransverseMercatorExact` are skipped.

To regenerate (after running `script/download-test-data.sh`):

```sh
./script/excerpt-tmcoords.py
```

### Full fixture

Use `script/download-test-data.sh` and the `test_full` feature flag to run against the original `TMcoords.dat`.  The test internally skips rows the 6th-order Krueger series cannot resolve (segments 7..13).
