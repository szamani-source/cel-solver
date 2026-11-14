#!/bin/bash
# WP 1.3a: 2D plane-strain elastic Lagrange-plus-remap prototype (pure Python, no OpenFOAM).
# Checks (prototype/lr2d.py --verify), material streaming at 0.25 c_p along the diagonal:
#   1. smooth P-waves along x, y and at 45 degrees, S-waves along x and at 45 degrees:
#      speed within 1 % of c_p or c_s, L2 error order >= 1.7 over 100/200/400 cells per metre
#   2. a P-wave along x matches prototype/lr1d.py at 200 cells to 1e-10
#   3. a problem mirrored about the diagonal gives mirrored fields to 1e-12
#   4. step P-wave along x: no new extrema above 1 % of the jump; mass and both momentum
#      components conserved to 1e-12. The 45 degree step is run and reported (decision 0012).
#   5. fixed-wall box keeps its mass
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
python3 "$root/prototype/lr2d.py" --verify
