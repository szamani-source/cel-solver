#!/bin/bash
# WP 3.4 E sweep (decision 0034): T0 tank discharge rate at E = 0.1, 1 and
# 10 MPa, within 10 percent of the 0.1 MPa rate. FAILS on the current code
# (10 mm: -11 % at 1 MPa; 5 mm: -29 %); it is the done-when of the follow-up
# package and is not in CI. About 8 minutes at 10 mm on one core.
#
#   tests/manual/test_discharge_stiffness.sh                       # 10 mm, serial
#   tests/manual/test_discharge_stiffness.sh --cell 0.005 --procs 4 0.1 1
#
# With rule 12 (WP 3.4b, decision 0035): prefix CEL_EXTRA_MATERIAL='normalCap yes;'
# (10 mm then passes: 78.6 / 73.8 / 84.8; 5 mm still fails: 82.6 / 69.7).
set -eu
root="$(cd "$(dirname "$0")/../.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/discharge/discharge_checks.py" "$work" --sweep "$@"
