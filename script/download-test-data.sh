#!/bin/sh

set -ex

REPO_ROOT=$(git rev-parse --show-toplevel)
OUTPUT_DIR="${REPO_ROOT}/test_fixtures/test_data_unzipped"

mkdir -p "${OUTPUT_DIR}"

curl -L https://sourceforge.net/projects/geographiclib/files/testdata/GeodTest.dat.gz | gunzip > "${OUTPUT_DIR}/GeodTest.dat"
curl -L https://sourceforge.net/projects/geographiclib/files/testdata/GeodTest-short.dat.gz | gunzip > "${OUTPUT_DIR}/GeodTest-short.dat"

