#!/bin/bash
# WP 2.1: Eulerian library skeleton (src/celEulerian, class eulerianSolid).
#
# A. smokeBox with no time directory: every field is created from defaults
#    (rho from materialProperties) and written; time step and mass match the
#    hand calculation.
# B. smokeBox with setFields (half the box rho 2000 and U (1 0 0)): mass and
#    momentum at the start match the hand calculation to 1e-12 relative.
# C. Case B decomposed on 2 processors prints the same time steps, mass and
#    momentum as the serial run (1e-12 relative).
#
# WP 2.3 (decision 0016): since the Lagrange step is in, case B moves (the moving
# half meets the fixed walls), so its time step, mass and momentum change during
# the run. Case A (at rest) still checks that every step and value stays constant;
# B checks the start values against the hand calculation; C now compares every
# printed mass and momentum with the serial run, not only the last.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

# Hand calculation: dt = maxCo * dx / c, c = sqrt((K + 4G/3)/rho_min)
dtHand=$(python3 -c "
E=1e7; nu=0.3; rho=1800
K=E/(3*(1-2*nu)); G=E/(2*(1+nu))
print(repr(0.5*0.1/((K+4*G/3)/rho)**0.5))")
echo "hand-calculated deltaT = $dtHand"

# check.py LOG MASS MX MY MZ [moving]: run ended, the stable time step at the
# start matches dtHand; at rest (default): every printed deltaT matches dtHand
# and every printed mass and momentum matches; moving: the first printed (start)
# mass and momentum match
check() {
    python3 - "$@" <<PYEOF
import re, sys
log, mass, mx, my, mz = sys.argv[1], *map(float, sys.argv[2:6])
moving = len(sys.argv) > 6 and sys.argv[6] == 'moving'
txt = open(log).read()
assert re.search(r'^End', txt, re.M), 'celFoam did not reach End'
dt = float(re.search(r'^Stable time step *deltaT = (\S+)', txt, re.M).group(1))
steps = [float(x) for x in re.findall(r'^Time = \S+  deltaT = (\S+)', txt, re.M)]
masses = [float(x) for x in re.findall(r'^Total mass *= (\S+)', txt, re.M)]
moms = re.findall(r'^Total momentum *= \((\S+) (\S+) (\S+)\)', txt, re.M)
assert steps and masses and moms, 'missing printed output'
hand = $dtHand
def rel(a, b, scale): return abs(a - b)/scale
print(f'{log}: deltaT {dt!r}, {len(steps)} steps, mass {masses[-1]!r}, momentum {moms[-1]}')
assert rel(dt, hand, hand) < 1e-12, f'deltaT {dt} != {hand}'
scale = max(abs(mx), abs(my), abs(mz), 1.0)
checked = [0] if moving else range(len(masses))
if not moving:
    # per-step deltaT is printed with 6 digits; it may only change by rounding
    assert all(rel(s, hand, hand) < 1e-5 for s in steps), 'per-step deltaT changed'
for i in checked:
    assert rel(masses[i], mass, mass) < 1e-12, f'mass {masses[i]} != {mass}'
    for got, want in zip(map(float, moms[i]), (mx, my, mz)):
        assert rel(got, want, scale) < 1e-12, f'momentum {moms[i]} != {(mx, my, mz)}'
PYEOF
}

# fields DIR: the latest time directory under DIR holds all five fields
fields() {
    last=$(cd "$1" && ls -d [0-9]* | sort -g | tail -n 1)
    for f in rho U sigma epsP epsPEq
    do
        [ -f "$1/$last/$f" ] || { echo "FAIL: $1/$last/$f not written"; exit 1; }
    done
    echo "$1/$last holds rho U sigma epsP epsPEq"
}

echo "== A. defaults, no time directory"
cp -r "$root/tutorials/smokeBox" "$work/A"
cd "$work/A"
blockMesh > log.blockMesh 2>&1
celFoam > log.celFoam 2>&1
check log.celFoam 1800 0 0 0
fields .
last=$(ls -d [0-9]* | sort -g | tail -n 1)
grep -q '^internalField *uniform 1800;' "$last/rho" \
    || { echo "FAIL: default rho not uniform 1800"; exit 1; }

echo "== B. setFields, serial"
cp -r "$root/tutorials/smokeBox" "$work/B"
cd "$work/B"
cp -r 0.orig 0
blockMesh > log.blockMesh 2>&1
setFields > log.setFields 2>&1
cp -r . "$work/C"
celFoam > log.celFoam 2>&1
check log.celFoam 1900 1000 0 0 moving
fields .

echo "== C. setFields, 2 processors"
cd "$work/C"
decomposePar > log.decomposePar 2>&1
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
mpirun -np 2 --oversubscribe celFoam -parallel > log.celFoam 2>&1 \
    || { tail -n 30 log.celFoam; exit 1; }
check log.celFoam 1900 1000 0 0 moving
fields processor0
fields processor1

# Serial and parallel print the same times and time steps, and the same
# number of writes (values were checked against the hand calculation above)
diff <(grep '^Time = ' "$work/B/log.celFoam") <(grep '^Time = ' "$work/C/log.celFoam") \
    || { echo "FAIL: time steps differ between serial and parallel"; exit 1; }
a=$(grep -c '^Total mass' "$work/B/log.celFoam")
b=$(grep -c '^Total mass' "$work/C/log.celFoam")
[ "$a" = "$b" ] || { echo "FAIL: $a writes serial vs $b parallel"; exit 1; }
python3 - "$work/B/log.celFoam" "$work/C/log.celFoam" <<'PYEOF'
import re, sys
def values(log):
    txt = open(log).read()
    m = [float(x) for x in re.findall(r'^Total mass *= (\S+)', txt, re.M)]
    p = [tuple(map(float, x)) for x in
         re.findall(r'^Total momentum *= \((\S+) (\S+) (\S+)\)', txt, re.M)]
    return m, p
(ms, ps), (mp, pp) = values(sys.argv[1]), values(sys.argv[2])
for a, b in zip(ms, mp):
    assert abs(a - b) <= 1e-12*abs(a), f'mass serial {a} parallel {b}'
scale = max(max(abs(c) for c in p) for p in ps)
for a, b in zip(ps, pp):
    assert max(abs(x - y) for x, y in zip(a, b)) <= 1e-12*scale, \
        f'momentum serial {a} parallel {b}'
print(f'{len(ms)} printed masses and momenta agree to 1e-12')
PYEOF
echo "serial and parallel agree"
