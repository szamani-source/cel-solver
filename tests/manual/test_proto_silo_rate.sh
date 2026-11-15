#!/bin/bash
# WP 3.3p done-when rate checks, open in WP 3.3q (decision 0024). Not in CI: fails now.
# prototype/discharge2d.py --rate-check (about 10 minutes): 10 mm cells, E = 0.1 MPa,
# D = 0.04 / 0.06 / 0.08 m and H = 0.9 m: steady rate (0.1 s windows within 10 %),
# exponent 1.5 +- 0.2, C in 1.0 to 1.9 (docs/spec/silo-reference.md section 6),
# corner speed below 1 % of the outlet speed, a crater, head independence 10 %.
set -eu
root="$(cd "$(dirname "$0")/../.." && pwd)"
python3 -W ignore "$root/prototype/discharge2d.py" --rate-check
