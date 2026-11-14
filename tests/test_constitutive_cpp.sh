#!/bin/bash
# WP 2.2: C++ constitutive library (src/celConstitutive). WP 2.5a added Drucker-Prager,
# WP 2.5b Mohr-Coulomb, WP 2.5c its tension cut-off (sigmaT).
#
# 1. Kernel: builds the plain C++ kernel with g++ (no OpenFOAM needed), unit-tests the
#    symmetric 3x3 eigen-solver (repeated eigenvalues included), replays the WP 1.4
#    Python reference paths (uniaxial strain load-unload, triaxial compression
#    and extension, simple shear, 1000 random increments) for linear elastic,
#    von Mises (perfect and hardening), Drucker-Prager (cohesive and cohesionless,
#    plus the extension and plane-strain fits) and Mohr-Coulomb (cohesive and
#    cohesionless, and with tension cut-offs sigma_t = 0 and 0.5 c cot(phi):
#    triaxial and hydrostatic tension, tension then shear, uniaxial strain, 1000
#    tension-biased random increments) and requires stress and history to match to a relative 1e-10, with
#    apex returns on the Drucker-Prager random paths and, on each Mohr-Coulomb path,
#    the edge or apex return it is there to test.
# 2. OpenFOAM wrapper: if OpenFOAM is loaded, checks libcelConstitutive.so was built by
#    Allwmake, builds and runs tests/constitutive/wrapperCheck (run-time selection,
#    symmTensor component order, von Mises and Drucker-Prager returns, the fit word)
#    and checks that celFoam reports the model selected in the smokeBox tutorial,
#    then again with druckerPrager selected, and with mohrCoulomb plus sigmaT
#    (and that an out-of-range sigmaT stops celFoam).
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

# WP 2.5a: the same case with the Drucker-Prager model selected
foamDictionary -entry type -set druckerPrager constant/materialProperties > /dev/null
foamDictionary -entry phi -add 30 constant/materialProperties > /dev/null
foamDictionary -entry c -add 10e3 constant/materialProperties > /dev/null
foamDictionary -entry psi -add 5 constant/materialProperties > /dev/null
foamDictionary -entry fit -add compression constant/materialProperties > /dev/null
celFoam > log.celFoam.dp 2>&1 || { tail -n 30 log.celFoam.dp; exit 1; }
grep '^Constitutive model' log.celFoam.dp
grep -q '^Constitutive model *druckerPrager$' log.celFoam.dp
grep -q '^End' log.celFoam.dp
echo "celFoam selected druckerPrager"

# WP 2.5c: Mohr-Coulomb with the optional tension cut-off sigmaT
foamDictionary -entry type -set mohrCoulomb constant/materialProperties > /dev/null
foamDictionary -entry fit -remove constant/materialProperties > /dev/null
foamDictionary -entry sigmaT -add 5e3 constant/materialProperties > /dev/null
celFoam > log.celFoam.mc 2>&1 || { tail -n 30 log.celFoam.mc; exit 1; }
grep '^Constitutive model' log.celFoam.mc
grep -q '^Constitutive model *mohrCoulomb$' log.celFoam.mc
grep -q '^End' log.celFoam.mc
echo "celFoam selected mohrCoulomb with sigmaT"

# An out-of-range sigmaT (above c cot(phi) = 17.3 kPa) must stop celFoam
foamDictionary -entry sigmaT -set 1e5 constant/materialProperties > /dev/null
if celFoam > log.celFoam.bad 2>&1
then
    echo "FAIL: celFoam accepted sigmaT > c cot(phi)"
    exit 1
fi
grep -q 'sigmaT must be in' log.celFoam.bad
echo "celFoam rejects sigmaT > c cot(phi)"
