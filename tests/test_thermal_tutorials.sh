#!/bin/bash
# WP T.7 / decision 0140: the solid-fluid-temperature tutorials
# (tutorials/solidFluidThermal).
#
# [1] Layout: tutorials/solidFluidThermal/README.md lists every tutorial folder,
#     and every folder has an executable Allrun and Allclean and a README.md with
#     the sections Setup, Numerics, How to run and Expected results; every figure
#     a README shows exists, and each README shows at least one.
# [2] Short runs from a copy (about 40 s in all): hotAirPackedBed to 0.5 s,
#     movingBedHeatedWall to 0.1 s, heatedAirChannel to 0.02 s end with "End",
#     no NaN, and the energy balance in log.celFoam under 1e-10 of the thermal
#     energy. hotSiloDischarge runs in tests/test_hot_silo_tutorial.sh.
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

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

echo "[2] short runs"
short() {    # tutorial, end time
    local t="$1" te="$2"
    cp -r "$tut/$t" "$work/$t"
    foamDictionary -entry endTime -set "$te" "$work/$t/system/controlDict" > /dev/null 2>&1
    foamDictionary -entry writeInterval -set "$te" "$work/$t/system/controlDict" > /dev/null 2>&1
    (cd "$work/$t" && ./Allrun > log.Allrun 2>&1)
    local log="$work/$t/log.celFoam"
    if [ ! -f "$log" ] || ! grep -q '^End' "$log"; then
        bad "$t: celFoam did not finish"; tail -n 20 "$log" 2>/dev/null; return
    fi
    grep -qi ' nan' "$log" && bad "$t: NaN in log.celFoam"
    python3 - "$log" "$t" "$te" <<'PY' || fail=1
import re, sys
log = open(sys.argv[1]).read()
e0 = float(re.search(r"^Thermal energy at start\s*=\s*(\S+)", log, re.M).group(1))
bal = max(abs(float(v)) for v in re.findall(r"^Energy balance\s*=\s*(\S+)", log, re.M))
last = float(re.findall(r"^Time = (\S+)", log, re.M)[-1])
ok = bal <= 1e-10 * e0 and last >= 0.95 * float(sys.argv[3])
print(f"  {sys.argv[2]}: ran to t = {last:g}, energy balance {bal / e0:.1e} of the thermal energy"
      + ("" if ok else ": FAIL"))
sys.exit(0 if ok else 1)
PY
}
short hotAirPackedBed 0.5
short movingBedHeatedWall 0.1
short heatedAirChannel 0.02

if [ "$fail" -eq 0 ]; then echo "test_thermal_tutorials: PASS"; else echo "test_thermal_tutorials: FAIL"; fi
exit "$fail"
