#!/bin/bash
# WP 0.1: celFoam runs on a 10^3 hex box and uses the expected stable time step.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
cp -r "$root/tutorials/solid/smokeBox" "$work/"
cd "$work/smokeBox"

blockMesh > log.blockMesh 2>&1
celFoam > log.celFoam 2>&1

grep -q '^End' log.celFoam

# Expected deltaT = maxCo * dx / c, c = sqrt((K + 4G/3)/rho)
expected=$(python3 -c "
E=1e7; nu=0.3; rho=1800
K=E/(3*(1-2*nu)); G=E/(2*(1+nu))
c=((K+4*G/3)/rho)**0.5
print(0.5*0.1/c)")
actual=$(sed -n 's/^Stable time step *deltaT = //p' log.celFoam)
python3 -c "
import sys
e=$expected; a=$actual
print('expected', e, 'actual', a)
sys.exit(0 if abs(a-e)/e < 1e-6 else 1)"
