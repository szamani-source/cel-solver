#!/bin/bash
# WP 3.2s done-when (column collapse, from WP 3.2q): 2D granular column collapse,
# a = 0.5 and 2, T0 Mohr-Coulomb (c = 0, phi 35), tolerances from
# docs/spec/silo-reference.md section 6 (decision 0019): a = 0.5 keeps its top (within
# 1 cell); run-out within 0.5 to 1.5 times 1.2 a; deposit slope at most 37 deg as a
# line fitted from 10 to 90 % of the final height (decision 0025, chosen by Sia
# 2026-11-15; single-cell slopes are printed for information); mass to 1e-10.
# Takes about 4 minutes.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
python3 -W ignore "$root/prototype/void2d.py" --collapse
