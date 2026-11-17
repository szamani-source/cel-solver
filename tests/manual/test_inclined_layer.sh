#!/bin/bash
# WP 3.6b done-when (decision 0039): steady inclined layer (Bagnold profile)
# with mu(I): surface speed and profile within 5 % of the exact solution at
# E = 0.1 MPa, 5 mm cells; E = 1 and 10 MPa, 2.5 mm, and a smooth lid instead
# of the free surface recorded. FAILS on the current code (the flow is 5 to
# 7 times too slow at 0.1 MPa: numerical viscosity about 0.06 rho c_s h).
# It is the profile check of WP 3.6c. About 3 minutes on one core.
#
#   tests/manual/test_inclined_layer.sh               # all cases
#   tests/manual/test_inclined_layer.sh 0.1:5 1:5:lid # E_MPa:cell_mm[:lid]
set -eu
root="$(cd "$(dirname "$0")/../.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/layer/inclined_layer.py" "$work" "$@"
