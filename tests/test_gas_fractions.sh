#!/bin/bash
# WP F.1: gas fractions and the celFluid skeleton (decisions 0125, 0126).
# A settled 2D column with a void top runs in celFoam with and without
# constant/gasProperties; tests/gas/gas_fraction_checks.py checks phi_s against
# the set packing, k against tools/gas_reference.py, the solid fields bit for
# bit against the runs without gas, and 2 processors against serial.
# About 1 minute.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

python3 "$root/tests/gas/gas_fraction_checks.py" "$work"
