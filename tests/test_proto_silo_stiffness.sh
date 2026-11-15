#!/bin/bash
# WP 3.3r done-when: the T0 discharge rate does not depend on the stiffness
# (decision 0028). prototype/discharge2d.py --stiffness-check: 10 mm cells, D = 0.06 m,
# H = 0.6 m, E = 0.1, 1 and 10 MPa (T0) in three parallel processes. Checks: rates
# within 10 % of each other, 0.1 s windows within 10 %, stagnant corners (mean speed
# over 0.3 s) below 1 % of the outlet speed, nothing through the closed outlet,
# mass balance 1e-8. Before rules 12 and 13 of void2d.py: 81 / 73 / 54 kg/s/m.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
python3 -W ignore "$root/prototype/discharge2d.py" --stiffness-check
