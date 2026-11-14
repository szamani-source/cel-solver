#!/bin/bash
# WP 2.4b: the full C++ cycle (eulerianSolid::lagrangeStep then remap, as celFoam
# runs them) reproduces the 2D prototype checks (decisions 0013, 0017, 0018, 0022).
#
# Builds the test driver tests/eulerian/cycleCheck (celCycleCheck: initial U and
# sigma from 0/, fixed dt, records totals and wave-frame ranges every step) and runs
# tests/eulerian/cycle_checks.py:
#
# - Elastic P and S waves along x and at 45 degrees, material streaming diagonally
#   at 0.25 c_p: speed within 1 %, order >= 1.7 (25, 50, 100 cells), step waves no
#   new extrema above 1 % after 10 cells of travel, mass and momentum to 1e-12.
# - Pre-stressed disc in prescribed rigid rotation, one turn, 100 cells across:
#   principal stresses within 2 %, below 1 % beyond the 6-cell smeared edge.
# - von Mises step loading on a 1000 x 2 x 2 3D mesh, streaming: precursor and
#   plastic wave speeds and the HEL within 2 % of the analytic values.
# - 45 degree impact (von Mises, Drucker-Prager, elastic): beyond the along-x range
#   by at most 1 % after the start-up; 2 processors match serial to 1e-12.
# - Stable at the default maxCo (read from eulerianSolid.C), 2D and 3D, 4 crossings.
#
# Takes about 3 minutes. `cycle_checks.py WORKDIR --study` scans maxCo (decision 0022).
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/eulerian/cycleCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

maxCo=$(sed -n 's/.*getOrDefault<scalar>("maxCo", \([0-9.]*\)).*/\1/p' \
    "$root/src/celEulerian/eulerianSolid/eulerianSolid.C")
[ -n "$maxCo" ] || { echo "FAIL: default maxCo not found in eulerianSolid.C"; exit 1; }

python3 "$root/tests/eulerian/cycle_checks.py" "$work" "$maxCo"
