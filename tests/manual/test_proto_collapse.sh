#!/bin/bash
# WP 3.2q done-when (split from WP 3.2p, decision 0021): 2D granular column collapse,
# a = 0.5 and 2, T0 Mohr-Coulomb (c = 0, phi 35), tolerances from
# docs/spec/silo-reference.md section 6 (decision 0019): a = 0.5 keeps its top (within
# 1 cell); run-out within 0.5 to 1.5 times 1.2 a; no deposit slope steeper than 37 deg.
# Not yet in tests/run_tests.sh: run-outs pass since WP 3.2q (decision 0023), the deposit
# slope still fails (WP 3.2r). Move it to tests/ when it passes. Takes about 4 minutes.
set -eu
root="$(cd "$(dirname "$0")/../.." && pwd)"
python3 -W ignore "$root/prototype/void2d.py" --collapse
