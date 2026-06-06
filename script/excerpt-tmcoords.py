#!/usr/bin/env python3

"""Sample Karney's TMcoords.dat to produce the built-in TMcoords excerpt fixture.

`TMcoords.dat` (287_000 rows; downloaded by `script/download-test-data.sh`)
is organised in 13 segments described at:
  https://geographiclib.sourceforge.io/C++/doc/transversemercator.html#testmerc
"""

import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SRC = REPO_ROOT / "test_fixtures" / "test_data_unzipped" / "TMcoords.dat"
OUT = REPO_ROOT / "test_fixtures" / "TMcoords-excerpt.dat"

# Restrict samples to the practical Gauss-Schreiber eta range, well within
# the Krueger 6th-order series' FP64-accurate domain.
ETA_MAX = 0.5

# Row spans (1-indexed inclusive) of the 13 segments in the original TMcoords.dat.
SEGMENTS = {
    1: (1, 250_000),
    2: (250_001, 251_000),
    3: (251_001, 252_000),
    4: (252_001, 253_000),
    5: (253_001, 254_000),
    6: (254_001, 255_000),
    # 7..13 require TransverseMercatorExact
    # 7: (255_001, 256_000),
    # 8: (256_001, 258_000),
    # 9:  (258_001, 283_000),
    # 10: (283_001, 284_000),
    # 11: (284_001, 285_000),
    # 12: (285_001, 286_000),
    # 13: (286_001, 287_000),
}

# 10 x 5 stratified lat x lon grid for segment 1.
MAIN_GRID = (10, 5)
# Rows sampled from each non-segment-1 entry in SEGMENTS.
SEGMENT_PICKS = 10


def eta(lat: float, lon: float) -> float:
    """Approximate Gauss-Schreiber eta for (lat, lon) on WGS84."""
    sphi, cphi = math.sin(math.radians(lat)), math.cos(math.radians(lat))
    slam, clam = math.sin(math.radians(lon)), math.cos(math.radians(lon))
    if cphi == 0.0:
        return 0.0
    h = math.hypot(sphi / cphi, clam)
    return math.asinh(abs(slam) / h) if h else math.inf


def main() -> int:
    if not SRC.exists():
        sys.stderr.write(
            f"error: {SRC} not found; run script/download-test-data.sh first\n"
        )
        return 1

    main_grid: dict[tuple[int, int], str] = {}
    seg_rows: dict[int, list[str]] = {seg: [] for seg in SEGMENTS if seg != 1}

    with SRC.open() as f:
        for line_idx, line in enumerate(f, 1):
            line = line.rstrip("\n")
            lat, lon = (float(x) for x in line.split()[:2])
            if eta(lat, lon) >= ETA_MAX:
                continue
            for seg, (start, end) in SEGMENTS.items():
                if start <= line_idx <= end:
                    if seg == 1:
                        i = min(int(lat / (90 / MAIN_GRID[0])), MAIN_GRID[0] - 1)
                        j = min(int(lon / (90 / MAIN_GRID[1])), MAIN_GRID[1] - 1)
                        main_grid.setdefault((i, j), line)
                    else:
                        seg_rows[seg].append(line)
                    break

    rows = [
        main_grid[(i, j)]
        for i in range(MAIN_GRID[0])
        for j in range(MAIN_GRID[1])
        if (i, j) in main_grid
    ]
    for buf in seg_rows.values():
        n = min(SEGMENT_PICKS, len(buf))
        rows.extend(buf[i * len(buf) // n] for i in range(n))

    OUT.write_text("\n".join(rows) + "\n")
    print(f"Wrote {len(rows)} rows to {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
