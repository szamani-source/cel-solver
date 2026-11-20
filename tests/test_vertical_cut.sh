#!/bin/bash
# WP C.2: a cohesive column with a free vertical face stands or fails at the
# critical cohesion of a vertical cut (Terzaghi), the same on 5 and 2.5 mm
# cells, Mohr-Coulomb and Drucker-Prager (decision 0112). About 4 minutes.
#
# Builds the test driver celVoidCheck and runs tests/cohesion/vertical_cut.py.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT
mkdir -p "$work"

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/cohesion/vertical_cut.py" "$work" "$@"
