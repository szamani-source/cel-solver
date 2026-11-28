#!/bin/bash
# WP TU.1: the solid-only tutorials (tutorials/solid, decision 0140).
#
# [1] Layout: tutorials/solid/README.md lists every tutorial folder, and every
#     folder has an executable Allrun and Allclean and a README.md with the
#     sections Setup, Numerics, How to run and Expected results; every figure a
#     README shows exists, and each README shows at least one.
# [2] Each tutorial starts and runs a short stretch from a copy (soft material,
#     short end time; about 3 minutes in all): the OpenFOAM-case tutorials end
#     with "End" past the outlet opening and no NaN; coutoHopper builds its
#     7,584-cell mesh (Mesh OK); setupPlaneTank writes, checks and runs its
#     case to 0.02 s. stripFooting runs the case of tests/test_footing.sh, so it
#     is not run again here.
set -u
root="$(cd "$(dirname "$0")/.." && pwd)"
tut="$root/tutorials/solid"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
fail=0
bad() { echo "FAIL: $*"; fail=1; }

echo "[1] layout and READMEs"
[ -f "$tut/README.md" ] || bad "tutorials/solid/README.md missing"
for d in "$tut"/*/; do
    t="$(basename "$d")"
    grep -q "\`$t\`" "$tut/README.md" || bad "$t not listed in tutorials/solid/README.md"
    for f in Allrun Allclean; do
        [ -x "$d/$f" ] || bad "$t/$f missing or not executable"
    done
    r="$d/README.md"
    if [ ! -f "$r" ]; then bad "$t/README.md missing"; continue; fi
    for s in "## Setup" "## Numerics" "## How to run" "## Expected results"; do
        grep -q "^$s" "$r" || bad "$t/README.md has no '$s' section"
    done
    n=0
    for img in $(grep -o '!\[[^]]*\]([^)]*)' "$r" | sed 's/.*(\(.*\))/\1/'); do
        n=$((n + 1))
        [ -f "$d/$img" ] || bad "$t/README.md shows $img, which is missing"
    done
    [ "$n" -ge 1 ] || bad "$t/README.md shows no figure"
    echo "  $t: README sections and $n figure(s) checked"
done

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }
wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

echo "[2] short runs"
ok_log() {   # case dir, label, end time it must pass
    local log="$1/log.celFoam"
    if [ ! -f "$log" ] || ! grep -q '^End' "$log"; then
        bad "$2: celFoam did not finish"; tail -n 20 "$log" 2>/dev/null; return
    fi
    grep -qi ' nan' "$log" && bad "$2: NaN in log.celFoam"
    last=$(grep '^Time = ' "$log" | tail -n 1 | awk '{print $3}')
    awk -v a="$last" -v b="$3" 'BEGIN { exit !(a >= 0.95 * b) }' \
        || bad "$2: ended at t = $last, before $3"
    echo "  $2: ran to t = $last"
}
short() {    # tutorial, end time, Allrun arguments...
    local t="$1" te="$2"; shift 2
    cp -r "$tut/$t" "$work/$t"
    foamDictionary -entry endTime -set "$te" "$work/$t/system/controlDict" > /dev/null 2>&1
    (cd "$work/$t" && ./Allrun "$@" > log.Allrun 2>&1)
}

short smokeBox 0.01 && ok_log "$work/smokeBox" smokeBox 0.01
short rigidBodyMotion 0.01 && ok_log "$work/rigidBodyMotion" rigidBodyMotion 0.01
for t in siloDischarge siloFed stlSilo stlHopper; do
    short $t 0.35 0.1 && ok_log "$work/$t" $t 0.35
done
cp -r "$tut/siloCohesion" "$work/siloCohesion"
foamDictionary -entry endTime -set 0.35 "$work/siloCohesion/system/controlDict" > /dev/null 2>&1
foamDictionary -entry E -set 1e6 "$work/siloCohesion/constant/materialProperties" > /dev/null 2>&1
(cd "$work/siloCohesion" && ./Allrun 140 > log.Allrun 2>&1)
ok_log "$work/siloCohesion/D140" siloCohesion 0.35

cp -r "$tut/coutoHopper" "$work/coutoHopper"
(cd "$work/coutoHopper" && CELFOAM_DIR="$root" ./Allrun mc mesh > mesh.log 2>&1)
if grep -q "Mesh OK" "$work/coutoHopper/mesh.log" && grep -q "cells: *7584" "$work/coutoHopper/mesh.log"
then echo "  coutoHopper: 7,584-cell mesh, Mesh OK"
else bad "coutoHopper mesh"; tail -n 20 "$work/coutoHopper/mesh.log"
fi

cp -r "$tut/setupPlaneTank" "$work/setupPlaneTank"
(cd "$work/setupPlaneTank" && sed -i "s#^tools=.*#tools=$root/tools#" Allrun \
    && ./Allrun 0.02 > log.Allrun 2>&1)
ok_log "$work/setupPlaneTank/case" setupPlaneTank 0.02

if [ "$fail" -eq 0 ]; then echo "solid tutorials: PASS"; else echo "solid tutorials: FAILED"; fi
exit $fail
