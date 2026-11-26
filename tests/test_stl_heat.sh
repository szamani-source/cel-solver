#!/bin/bash
# WP T.2b: temperature with immersed STL walls (decision 0135).
# tests/thermal/stl_heat_checks.py: tutorials/stlHopper with a hot lower
# layer: energy balance 1e-12 with STL walls and an outlet, material within
# the start's range, uniform temperature kept to 1e-12, heat carried out
# c T times the discharged mass. About 80 seconds.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

python3 "$root/tests/thermal/stl_heat_checks.py" "$work"
