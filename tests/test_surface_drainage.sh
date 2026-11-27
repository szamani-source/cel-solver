#!/bin/bash
# WP F.8i: the pore pressure drains at the material's own surface
# (gasProperties frictionOnly; decision 0133). tests/gas/surface_drainage_checks.py:
# a still layer with a partly filled top cell against the series solution,
# 2 processors against serial.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

python3 "$root/tests/gas/surface_drainage_checks.py" "$work"
