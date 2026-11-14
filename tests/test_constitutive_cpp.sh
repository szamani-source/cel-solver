#!/bin/bash
# WP 2.2 / 2.5b: C++ constitutive library (src/celConstitutive).
#
# 1. Kernel: builds the plain C++ kernel with g++ (no OpenFOAM needed), unit-tests the
#    symmetric 3x3 eigen-solver (repeated eigenvalues included), replays the WP 1.4
#    Python reference paths (uniaxial strain load-unload, triaxial compression
#    and extension, simple shear, 1000 random increments) for linear elastic,
#    von Mises (perfect and hardening) and Mohr-Coulomb (cohesive and cohesionless;
#    each path must hit its edge or apex return) and requires stress and history to
#    match to a relative 1e-10.
# 2. OpenFOAM wrapper: if OpenFOAM is loaded, checks libcelConstitutive.so was built by
#    Allwmake, builds and runs tests/constitutive/wrapperCheck (run-time selection,
#    symmTensor component order, a von Mises return) and checks that celFoam reports the model selected in the
#    smokeBox tutorial.
#    Skipped (and said so) without OpenFOAM, except in CI, where it is required.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

kernel="$root/src/celConstitutive/kernel"

echo "== Kernel test (g++)"
g++ -std=c++17 -O2 -Wall -Wextra -Werror -I"$kernel" \
    "$root/tests/constitutive/compareReference.C" "$kernel"/*.C \
    -o "$work/compareReference"
g++ -std=c++17 -O2 -Wall -Wextra -Werror -I"$kernel" \
    "$root/tests/constitutive/eigenCheck.C" "$kernel/celEigen.C" \
    -o "$work/eigenCheck"
"$work/eigenCheck"
echo
python3 "$root/tests/gen_constitutive_reference.py" "$work/reference.txt"
"$work/compareReference" "$work/reference.txt"

echo
echo "== OpenFOAM wrapper"
if [ -z "${WM_PROJECT_DIR:-}" ]
then
    if [ -n "${CI:-}" ]
    then
        echo "FAIL: OpenFOAM environment not loaded in CI"
        exit 1
    fi
    echo "SKIPPED: OpenFOAM not loaded, library and celFoam not checked here"
    exit 0
fi

lib="$FOAM_USER_LIBBIN/libcelConstitutive.so"
[ -f "$lib" ] || { echo "FAIL: $lib not found (run ./Allwmake)"; exit 1; }
echo "found $lib"

wmake "$root/tests/constitutive/wrapperCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }
celConstitutiveWrapperCheck

cp -r "$root/tutorials/smokeBox" "$work/"
cd "$work/smokeBox"
blockMesh > log.blockMesh 2>&1
celFoam > log.celFoam 2>&1
grep '^Constitutive model' log.celFoam
grep -q '^Constitutive model *linearElastic$' log.celFoam
grep -q '^End' log.celFoam
echo "celFoam selected the smokeBox model"
