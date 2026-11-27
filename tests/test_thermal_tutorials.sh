#!/bin/bash
# WP T.7 / decision 0140: the solid-fluid-temperature tutorials
# (tutorials/solidFluidThermal).
#
# [1] Layout: tutorials/solidFluidThermal/README.md lists every tutorial folder,
#     and every folder has an executable Allrun and Allclean and a README.md with
#     the sections Setup, Numerics, How to run and Expected results; every figure
#     a README shows exists, and each README shows at least one.
# The short runs are separate: tests/test_hot_silo_tutorial.sh (hotSiloDischarge).
set -u
root="$(cd "$(dirname "$0")/.." && pwd)"
tut="$root/tutorials/solidFluidThermal"
fail=0
bad() { echo "FAIL: $*"; fail=1; }

echo "[1] layout and READMEs"
[ -f "$tut/README.md" ] || bad "tutorials/solidFluidThermal/README.md missing"
n_tut=0
for d in "$tut"/*/; do
    t="$(basename "$d")"
    n_tut=$((n_tut + 1))
    grep -q "\`$t\`" "$tut/README.md" || bad "$t not listed in tutorials/solidFluidThermal/README.md"
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
[ "$n_tut" -ge 1 ] || bad "no tutorial in tutorials/solidFluidThermal"

if [ "$fail" -eq 0 ]; then echo "test_thermal_tutorials: PASS"; else echo "test_thermal_tutorials: FAIL"; fi
exit "$fail"
