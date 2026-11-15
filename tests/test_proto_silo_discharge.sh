#!/bin/bash
# WP 3.3p: target T0 in the 2D prototype, settle with the outlet closed, then discharge
# (prototype/discharge2d.py --check; pure Python, about 10 minutes). Coarse version:
# 10 mm cells, E = 1 MPa, rate from 97 to 88 percent of the mass. Checks
# (docs/spec/silo-reference.md section 6, decisions 0019 and 0024):
#   1. nothing leaves before the outlet opens; mass in tank + discharged constant to 1e-8
#   2. steady rate: 0.1 s window rates within 10 % of the fitted rate
#   3. rate law over D = 0.04, 0.06, 0.08 m: exponent 1.3 to 1.7 (C reported only;
#      the C band is checked in the full 5 mm runs at T0's E, docs/results/silo-prototype.md)
#   4. head independence: H = 0.9 m and 0.6 m rates within 10 %
#   5. flow pattern (D = 0.06 m): corner speed below 1 % of the outlet speed, a crater
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
python3 -W ignore "$root/prototype/discharge2d.py" --check
