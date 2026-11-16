#!/bin/bash
# WP 3.6a: mu(I) friction (Jop law, optional Barker-Gray regularisation) in the
# Drucker-Prager and Mohr-Coulomb models, Python then C++ (decision 0037).
#
# Done when (docs/work-packages.md):
# 1. Python (prototype/constitutive/check_mu_i.py, about 30 s): steady simple shear at
#    fixed pressure gives mu(I) to a relative 1e-6 for I = 1e-4 to 1 (placeholder sand
#    with the plain Jop law: Drucker-Prager psi 0 and 10 deg, Mohr-Coulomb psi 0;
#    Jop's glass beads with the Barker-Gray branch: both models), and mu2 = tan(phi)
#    reproduces the rate-independent paths bit for bit.
# 2. C++ kernel: matches Python on 30 mu(I) paths (simple shear at fixed pressure,
#    triaxial compression, pull to the apex, 300 random increments; six parameter
#    sets) to a relative 1e-10, and with mu2 = tan(phi) gives the rate-independent C++
#    paths of tests/gen_constitutive_reference.py bit for bit.
# 3. OpenFOAM: celFoam reads the mu(I) entries and reports the law; Barker-Gray
#    regularisation is refused for mu_s = tan 35 deg (no well-posed range).
#    Skipped (and said so) without OpenFOAM, except in CI, where it is required.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
kernel="$root/src/celConstitutive/kernel"

echo "== 1. Python reference"
(cd "$root/prototype/constitutive" && python3 check_mu_i.py)

echo
echo "== 2. C++ kernel against Python"
g++ -std=c++17 -O2 -Wall -Wextra -Werror -I"$kernel" \
    "$root/tests/constitutive/compareReference.C" "$kernel"/*.C \
    -o "$work/compareReference"
python3 "$root/tests/gen_mu_i_reference.py" "$work/muI.txt"
"$work/compareReference" "$work/muI.txt"
echo
echo "-- mu2 = tan(phi): rate-independent C++ paths bit for bit"
python3 "$root/tests/gen_constitutive_reference.py" "$work/reference.txt"
"$work/compareReference" --mu2-equal-mus "$work/reference.txt" | grep -v '^ *returns'

echo
echo "== 3. OpenFOAM: celFoam with mu(I)"
if [ -z "${WM_PROJECT_DIR:-}" ]
then
    if [ -n "${CI:-}" ]
    then
        echo "FAIL: OpenFOAM environment not loaded in CI"
        exit 1
    fi
    echo "SKIPPED: OpenFOAM not loaded, celFoam not checked here"
    exit 0
fi

cp -r "$root/tutorials/smokeBox" "$work/"
cd "$work/smokeBox"
blockMesh > log.blockMesh 2>&1
mp=constant/materialProperties
set_entry() { foamDictionary -entry "$1" -add "$2" $mp > /dev/null 2>&1; }

# Jop's glass beads with the Barker-Gray branch: no elastic range, so every
# stressed cell goes through the mu(I) return
foamDictionary -entry type -set druckerPrager $mp > /dev/null 2>&1
set_entry phi 20.9
set_entry c 0
set_entry psi 0
set_entry mu2 0.6434
set_entry I0 0.279
set_entry d 0.53e-3
set_entry rhoS 2450
set_entry muIAlpha 1.9
celFoam > log.celFoam.bg 2>&1 || { tail -n 30 log.celFoam.bg; exit 1; }
grep '^Friction mu(I)' log.celFoam.bg
grep -q '^Friction mu(I).*Barker-Gray alpha 1.9, I_1^N 0.019' log.celFoam.bg
grep -q '^End' log.celFoam.bg
echo "celFoam ran druckerPrager with Barker-Gray mu(I)"

# Placeholder sand (decision 0037) with Mohr-Coulomb: plain Jop law
foamDictionary -entry type -set mohrCoulomb $mp > /dev/null 2>&1
foamDictionary -entry phi -set 35 $mp > /dev/null 2>&1
foamDictionary -entry mu2 -set 0.961813 $mp > /dev/null 2>&1
foamDictionary -entry d -set 1e-3 $mp > /dev/null 2>&1
foamDictionary -entry rhoS -set 2650 $mp > /dev/null 2>&1
foamDictionary -entry muIAlpha -remove $mp > /dev/null 2>&1
celFoam > log.celFoam.mc 2>&1 || { tail -n 30 log.celFoam.mc; exit 1; }
grep '^Friction mu(I)' log.celFoam.mc
grep -q '^Friction mu(I).*no regularisation' log.celFoam.mc
grep -q '^End' log.celFoam.mc
echo "celFoam ran mohrCoulomb with mu(I)"

# Barker-Gray needs a well-posed range of I; tan 35 deg has none
set_entry muIAlpha 1.9
if celFoam > log.celFoam.bad 2>&1
then
    echo "FAIL: celFoam accepted Barker-Gray regularisation for mu_s = tan 35 deg"
    exit 1
fi
grep -q 'no well-posed range' log.celFoam.bad
echo "celFoam refuses Barker-Gray for mu_s = tan 35 deg"

echo
echo "PASS"
