#!/bin/bash
# WP 3.6f (decision 0046): scatter band of the T0 discharge rate with mu(I)
# and the flow impedance. Runs each E three times with maxCo 0.40, 0.38,
# 0.36 (same physics, different round-off) and prints the rates, mean and
# range. A measurement, not a pass/fail check; compare bands, not single
# runs. About 8 minutes at 10 mm for 0.1 / 1 / 10 MPa on 3 cores.
#
#   tests/manual/test_discharge_scatter.sh                 # 10 mm, 0.1 1 10 MPa
#   tests/manual/test_discharge_scatter.sh --runs 5 1      # 5 runs at 1 MPa
set -eu
root="$(cd "$(dirname "$0")/../.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/discharge/scatter.py" "$work" "$@"
