#!/bin/bash
# WP 1.1: 1D elastic Lagrange-plus-remap prototype (pure Python, no OpenFOAM).
# Checks, all with the material streaming through the fixed mesh at 0.25 c:
#   1. smooth Gaussian pulse: L2 stress error converges at order >= 1.8 (200/400/800 cells)
#   2. square pulse: no new extrema (1e-6 of the jump), front speed within 1 % of sqrt(M/rho)
#   3. mass and momentum conserved to 1e-12
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
python3 "$root/prototype/lr1d.py" --verify
