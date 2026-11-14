#!/bin/bash
# WP 1.2: 1D elastic-plastic Lagrange-plus-remap prototype (pure Python, no OpenFOAM).
# von Mises without hardening, uniaxial strain, material streaming at 0.25 c_e:
#   1. step loading above the HEL: precursor speed within 1 % of sqrt((K+4G/3)/rho),
#      plastic wave within 2 % of sqrt(K/rho), precursor amplitude within 2 % of
#      Y (K+4G/3)/(2G)
#   2. no new extrema above 1 % of the jump after the impact start-up (decision 0008)
#   3. mass and momentum conserved to 1e-12
#   4. pre-strained patch carried 60 and 120 cells: total plastic strain kept to 1e-12,
#      none made by the model, edges 4 cells or less (10-90 %) and not widening
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
python3 "$root/prototype/lr1d_plastic.py" --verify
