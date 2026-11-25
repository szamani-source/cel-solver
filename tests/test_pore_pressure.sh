#!/bin/bash
# WP F.2: pore pressure in a rigid bed (decisions 0125, 0127).
# tests/gas/pore_pressure_checks.py: permeameter against Ergun at two flow
# rates and the Wen-Yu branch, step at a face against erfc at three times,
# gas mass in a sealed non-orthogonal box, 2 and 4 processors against
# serial, and the cost of gas on against off. About 2 minutes.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

python3 "$root/tests/gas/pore_pressure_checks.py" "$work"
