#!/bin/bash
# WP 3.6b/3.6c done-when (decisions 0039, 0040): steady inclined layer
# (Bagnold profile) with mu(I) and flowImpedance yes: surface speed and
# profile within 5 % of the exact solution at E = 0.1 MPa, 5 mm cells
# (checked); E = 1 and 10 MPa, 2.5 mm, and a smooth lid recorded. WP 3.6c2
# needs all three E within 5 % at 5 mm. FAILS now (0.1 MPa, 5 mm: 7.8 %,
# floor cells). --noflow runs the WP 3.6b code path. About 6 minutes.
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
