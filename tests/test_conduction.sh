#!/bin/bash
# WP T.1: conduction in the grains of a fixed bed (decision 0135).
# tests/thermal/conduction_checks.py: erfc and constant-flux faces against
# tools/thermal_reference.py (1 %), a linear temperature steady on cells leaning
# 27 deg (1e-10), energy conserved to 1e-12 with adiabatic walls, 2 and 4
# processors against serial (1e-10), temperature off leaves the solid byte for
# byte. About 1-2 minutes.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

python3 "$root/tests/thermal/conduction_checks.py" "$work"
