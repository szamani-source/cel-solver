#!/bin/bash
# WP 4.2: exact cells and faces covered by STL bodies (celBodyCells): a plane
# at 30 deg cell by cell, a faceted sphere at 10 and 20 cells per diameter
# (summed volume against the STL's own volume), a cube moving one cell and a
# turning cube, a box whose surface strips meet on a mesh face (decision
# 0093), and 2 processors against serial.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
python3 "$root/tests/bodyCells/cases.py" make "$work"
for c in plane sphere10 sphere20 cube cubeTurn seam; do
    (cd "$work/$c" && blockMesh > log.blockMesh 2>&1 \
        && celBodyCells > log.celBodyCells 2>&1) \
        || { tail -n 30 "$work/$c/log.celBodyCells"; exit 1; }
    grep '^Time' "$work/$c/log.celBodyCells" | sed "s/^/$c: /"
done
cd "$work/sphere20"
cp 0/alphaBody alphaBody.serial
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
decomposePar > log.decomposePar 2>&1
mpirun -np 2 --oversubscribe celBodyCells -parallel > log.parallel 2>&1 \
    || { tail -n 30 log.parallel; exit 1; }
reconstructPar -withZero -fields '(alphaBody)' > log.reconstruct 2>&1
mkdir -p parallel/0
mv 0/alphaBody parallel/0/alphaBody
mv alphaBody.serial 0/alphaBody
python3 "$root/tests/bodyCells/cases.py" check "$work"
