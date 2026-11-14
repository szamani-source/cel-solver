#!/bin/bash
# WP 3.1p: 2D tank under gravity, smooth and frictional walls, settling (pure Python).
# Checks (prototype/silo2d.py --verify), T0 material (Mohr-Coulomb c = 0, phi 35, psi 0):
#   1. vectorised Mohr-Coulomb (constitutive/batch.py) matches models.py to 1e-12 on
#      random paths, with main-plane, edge and apex returns all exercised
#   2. lake at rest: a column in exact discrete K0 equilibrium between smooth walls keeps
#      velocities below 1e-10 m/s for 1000 steps
#   3. smooth walls, settled from zero stress: bottom vertical stress within 2 % of
#      rho g H, s_xx / s_yy within 5 % of K0 = nu / (1 - nu)
#   4. frictional walls (20 deg), W = 0.4 m, H = 3 W: mean vertical stress within 10 %
#      of Janssen below one width from the top (K measured at the walls); wall shear /
#      normal within 5 % of tan(20 deg) where the wall slips; measured and fitted K
#      within 0.27 to 0.45 (decision 0019); weight = bottom force +
#      wall friction within 1 %
#   5. both settled cases: kinetic energy criterion reached, and the state unchanged over
#      a further 10 % of the run with damping off
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
python3 "$root/prototype/silo2d.py" --verify
