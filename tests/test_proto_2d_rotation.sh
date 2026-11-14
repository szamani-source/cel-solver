#!/bin/bash
# WP 1.3b: 2D plane-strain prototype with rotation and plasticity (pure Python, no OpenFOAM).
# Checks (prototype/lr2d_plastic.py --verify):
#   1. vectorised WP 1.4 models (constitutive/batch.py) match models.py to 1e-12
#   2. pre-stressed elastic disc in rigid rotation, one full turn, 100 cells across:
#      principal stresses and the stress tensor within 2 % of the rotated initial stress
#      inside the disc, principal direction within 1 degree, stress outside the disc
#      below 1 % (beyond 6 cells of smeared edge)
#   3. simple shear to shear strain 1: von Mises shear stress within 1 % of Y/sqrt(3) after
#      yield and never dropping (no oscillation); elastic Jaumann closed form within 1e-3 G;
#      Drucker-Prager under confinement stays on the cone
#   4. von Mises plastic step wave along x on the 2D mesh matches prototype/lr1d_plastic.py
#      (WP 1.2) at 1000 cells to 1e-9; mass and momentum to 1e-12
#   5. impact at 45 degrees (von Mises, Drucker-Prager, elastic) adds no extrema above 1 %
#      over the same impact along x, after the start-up (decision 0013)
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
python3 "$root/prototype/lr2d_plastic.py" --verify
