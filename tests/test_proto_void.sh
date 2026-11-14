#!/bin/bash
# WP 3.2p: one material plus void in the 2D prototype (free surface; pure Python).
# Checks (prototype/void2d.py --verify), decision 0021:
#   1. lake at rest with the free surface inside the domain (top cell half full): K0
#      column between smooth walls keeps velocities below 1e-8 m/s for 1000 steps,
#      no stress in empty cells
#   2. a block translating through void at constant velocity for 50 cells keeps its
#      mass to round-off, its edges within 3 cells and its stress at zero
#   3. a block falling freely under gravity accelerates at g within 0.5 % and lands
#      on a rigid floor without losing mass or passing through
#   4. a cohesionless block pulled apart separates (void opens inside the domain)
#      without carrying tension
#   5. 2D granular column collapse, a = 0.5 and 2, T0 Mohr-Coulomb: run-out within
#      0.5 to 1.5 times 1.2 a, a = 0.5 keeps its top (within 1 cell), no deposit
#      slope steeper than 37 deg (docs/spec/silo-reference.md section 6, decision 0019)
#   Mass conserved to 1e-10 in all cases.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
python3 -W ignore "$root/prototype/void2d.py" --verify
