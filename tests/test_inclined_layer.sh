#!/bin/bash
# WP 3.6b (decision 0039), CI part: an inclined layer under gravity not along
# a mesh axis stays at rest below the stopping angle (the face reconstruction
# fix), and the mu(I) layer at E = 0.1 MPa, 5 mm runs with mass kept (its
# speeds are checked against the WP 3.6c guard: surface within 10 %, mean
# speed at least 0.85 of exact, with flowImpedance yes). About 10 s.
#
# The WP 3.6b/3.6c done-when (Bagnold profile within 5 %) is
# tests/manual/test_inclined_layer.sh; it still fails at 5 mm (decision 0040).
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/layer/inclined_layer.py" "$work" --ci
