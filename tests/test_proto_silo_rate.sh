#!/bin/bash
# WP 3.3q done-when: the T0 discharge follows the 2D Beverloo law (decision 0027).
# prototype/discharge2d.py --rate-check: 5 mm cells, E = 0.1 MPa, D = 0.04 / 0.06 /
# 0.08 m at H = 0.6 m and D = 0.06 m at H = 0.9 m, four cases in parallel processes
# (about 20 to 25 minutes on four cores). Checks:
#   steady rate (0.1 s windows within 10 %), exponent 1.5 +- 0.2, C in 1.0 to 1.9
#   (docs/spec/silo-reference.md section 6), stagnant corners (mean speed over 0.3 s,
#   one elastic period of the bed) below 1 % of the outlet speed, a crater, head
#   independence 10 %, nothing through the closed outlet, mass balance 1e-8.
# At 10 mm the 0.04 m outlet is only 4 cells wide (C 0.75 against 0.98): hence 5 mm.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
python3 -W ignore "$root/prototype/discharge2d.py" --rate-check
