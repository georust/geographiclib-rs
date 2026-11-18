#!/bin/sh
#
# Regenerates `test_fixtures/geocentric_wgs84.dat` from the C++
# GeographicLib reference implementation.

set -ex

REPO_ROOT=$(git rev-parse --show-toplevel)

SRC="${REPO_ROOT}/script/generate-geocentric-fixture.cpp"
BIN=$(mktemp -t generate-geocentric-fixture.XXXXXX)
OUT="${REPO_ROOT}/test_fixtures/geocentric_wgs84.dat"

trap 'rm -f "${BIN}"' EXIT

EXTRA_FLAGS=""
if [ -n "${GEOGRAPHICLIB_PREFIX}" ]; then
    EXTRA_FLAGS="-I${GEOGRAPHICLIB_PREFIX}/include -L${GEOGRAPHICLIB_PREFIX}/lib"
fi

# shellcheck disable=SC2086
c++ -std=c++17 -O2 ${EXTRA_FLAGS} \
    "${SRC}" -lGeographicLib -o "${BIN}"

"${BIN}" > "${OUT}"
