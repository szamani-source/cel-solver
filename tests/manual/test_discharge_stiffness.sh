#!/bin/bash
# T0 tank discharge rate at E = 0.1, 1 and 10 MPa, within 10 percent of the
# 0.1 MPa rate. Done-when of WP 3.6d (decision 0045): with mu(I) and the flow
# impedance (--mui, the default here), at 10 mm and at 5 mm. FAILS on the
# current code: 10 mm 65.5 / 53.8 / 67.7 kg/s/m (1 MPa -18 %), 5 mm on two
# processors 77.6 / 61.2 (1 MPa -21 %). Not in CI. 10 mm takes about 5 minutes
# on one core (10 MPa is 4 of them); 5 mm 0.1 and 1 MPa about 10 minutes.
#
#   tests/manual/test_discharge_stiffness.sh                       # 10 mm, serial
#   tests/manual/test_discharge_stiffness.sh --cell 0.005 --procs 2 0.1 1
#   CEL_NO_MUI=1 tests/manual/test_discharge_stiffness.sh          # plain Mohr-Coulomb
#                                                                  # (WP 3.4, decision 0034)
# With rule 12 (WP 3.4b, decision 0035): prefix CEL_EXTRA_MATERIAL='normalCap yes;'.
set -eu
root="$(cd "$(dirname "$0")/../.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

mui="--mui"
[ -z "${CEL_NO_MUI:-}" ] || mui=""
python3 "$root/tests/discharge/discharge_checks.py" "$work" --sweep $mui "$@"
