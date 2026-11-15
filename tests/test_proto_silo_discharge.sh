#!/bin/bash
# WP 3.3p: target T0 in the 2D prototype, settle with the outlet closed, then discharge
# through the opened outlet (prototype/discharge2d.py --check; pure Python, about
# 3 minutes). Coarse: 10 mm cells, E = 0.1 MPa, D = 0.06 m, rate from 97 to 88 percent
# of the mass. Checks (decision 0026):
#   1. nothing leaves before the outlet opens
#   2. mass in tank + discharged mass constant to 1e-8
#   3. material flows out (positive rate and outlet speed)
#   4. head independence: rates at H = 0.9 m and 0.6 m within 10 %
# The rate law checks (steady rate, D^1.5, C, stagnant corners) fail for a numerical
# reason (decision 0026) and are in tests/manual/test_proto_silo_rate.sh for WP 3.3q.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
python3 -W ignore "$root/prototype/discharge2d.py" --check
