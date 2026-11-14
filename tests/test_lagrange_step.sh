#!/bin/bash
# WP 2.3: Lagrange step in C++ (eulerianSolid::lagrangeStep, decisions 0013, 0014, 0018).
#
# Builds the test driver tests/lagrange/lagrangeStepCheck (celLagrangeStepCheck: sets a
# state by hand and runs Lagrange steps only, no remap) and runs
# tests/lagrange/run_checks.py:
#
# - P- and S-waves (strain 1e-6) along x and at 45 degrees in a periodic 2D box: speed
#   within 1 %, smooth wave converges at order >= 1.7 over 25, 50, 100 cells; total
#   momentum conserved to 1e-12; rho V_L = m (decision 0014); first-order switch.
# - Stability at the default maxCo (read from eulerianSolid.C): 45 degree step in 2D,
#   body-diagonal step in 3D, no growth over 2 box crossings.
# - Prescribed rigid rotation of a uniform pre-stress, one turn: principal stresses
#   within 0.1 %, direction within 0.1 degree.
# - Prescribed simple shear to gamma 1: von Mises at Y/sqrt(3) within 0.1 %, no
#   oscillation; Drucker-Prager on its cone.
# - Step pulses reflecting from a fixed wall, a free surface and a symmetry plane:
#   amplitude and sign within 2 %.
# - 45 degree step wave on 2 and 4 processors matches serial to 1e-12.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/lagrange/lagrangeStepCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

maxCo=$(sed -n 's/.*getOrDefault<scalar>("maxCo", \([0-9.]*\)).*/\1/p' \
    "$root/src/celEulerian/eulerianSolid/eulerianSolid.C")
[ -n "$maxCo" ] || { echo "FAIL: default maxCo not found in eulerianSolid.C"; exit 1; }
echo "default maxCo = $maxCo"

python3 "$root/tests/lagrange/run_checks.py" "$work" "$maxCo"
