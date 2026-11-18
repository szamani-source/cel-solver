#!/bin/bash
# WP 4.1: rigid bodies read from STL move on their prescribed paths.
# tutorials/rigidBodyMotion: a fixed 45 degree wall, a cube pushed down with
# a velocity ramp (table), a cube translating and spinning about an edge.
# Checks (tests/rigidBody/check_motion.py): centre, rotation, velocity and
# omega against the exact path to 1e-12, the written STL to 1e-7, the fixed
# wall unchanged; 2 processors write the same poses as serial.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

cp -r "$root/tutorials/rigidBodyMotion" "$work/serial"
cd "$work/serial"
cp -r 0.orig 0
blockMesh > log.blockMesh 2>&1
celFoam > log.celFoam 2>&1
grep -q '^End' log.celFoam

cp -r "$root/tutorials/rigidBodyMotion" "$work/parallel"
cd "$work/parallel"
cp -r 0.orig 0
blockMesh > log.blockMesh 2>&1
decomposePar > log.decomposePar 2>&1
mpirun --allow-run-as-root --oversubscribe -np 2 celFoam -parallel \
    > log.celFoam 2>&1
grep -q '^End' log.celFoam

python3 "$root/tests/rigidBody/check_motion.py" "$work/serial" "$work/parallel"
