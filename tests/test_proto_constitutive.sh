#!/bin/bash
# WP 1.4: single-point constitutive driver (elastic, J2, Drucker-Prager, Mohr-Coulomb).
# Checks triaxial strengths against Mohr-Coulomb closed form, the Drucker-Prager
# compression fit, yield-surface consistency for 10,000 random increments (with apex
# and corner cases) and dilatancy in simple shear. Pure Python, no OpenFOAM needed.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
python3 "$root/prototype/constitutive/check_models.py"
