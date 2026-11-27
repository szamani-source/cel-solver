#!/bin/bash
# WP F.9: the gas tutorials (tutorials/solidFluid, decisions 0140, 0141).
#
# [1] Layout: tutorials/solidFluid/README.md lists every tutorial folder, and
#     every folder has an executable Allrun and Allclean and a README.md with
#     the sections Setup, Numerics, How to run and Expected results; every
#     figure a README shows exists, and each README shows at least one.
# [2] Short runs from a copy: airImpededHopper (with gas and dry/) to 0.2 s,
#     past the outlet opening; fluidisedDamBreak (fluidised, 10 mm) to 0.04 s,
#     whose fields equal those of tests/gas/collapse_cases.py's case to the
#     last digit. Each ends with "End" and no NaN. About 2 minutes.
set -u
root="$(cd "$(dirname "$0")/.." && pwd)"
tut="$root/tutorials/solidFluid"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
fail=0
bad() { echo "FAIL: $*"; fail=1; }

echo "[1] layout and READMEs"
[ -f "$tut/README.md" ] || bad "tutorials/solidFluid/README.md missing"
for d in "$tut"/*/; do
    t="$(basename "$d")"
    grep -q "\`$t\`" "$tut/README.md" || bad "$t not listed in tutorials/solidFluid/README.md"
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

echo "[2] short runs"
ok_log() {   # log, label, end time it must pass
    local log="$1"
    if [ ! -f "$log" ] || ! grep -q '^End' "$log"; then
        bad "$2: celFoam did not finish"; tail -n 20 "$log" 2>/dev/null; return
    fi
    grep -qi ' nan' "$log" && bad "$2: NaN in log.celFoam"
    last=$(grep '^Time = ' "$log" | tail -n 1 | awk '{print $3}')
    awk -v a="$last" -v b="$3" 'BEGIN { exit !(a >= 0.95 * b) }' \
        || bad "$2: ended at t = $last, before $3"
    echo "  $2: ran to t = $last"
}

cp -r "$tut/airImpededHopper" "$work/hopper"
sed -i 's/^endTime .*/endTime         0.2;/' "$work/hopper/system/controlDict"
(cd "$work/hopper" && ./Allrun > log.Allrun 2>&1)
ok_log "$work/hopper/log.celFoam" "airImpededHopper (gas)" 0.2
ok_log "$work/hopper/dry/log.celFoam" "airImpededHopper (dry)" 0.2
grep -q "^Pore gas: [0-9]* solves" "$work/hopper/log.celFoam" || bad "airImpededHopper: no gas solve"
grep -q "^Pore gas" "$work/hopper/dry/log.celFoam" && bad "airImpededHopper dry/: gas solved"

cp -r "$tut/fluidisedDamBreak" "$work/dam"
sed -i 's/^endTime .*/endTime         0.04;/' "$work/dam/system/controlDict"
(cd "$work/dam" && ./Allrun 10 1 gas > log.Allrun 2>&1)
ok_log "$work/dam/log.celFoam" "fluidisedDamBreak" 0.04
python3 - "$root" "$work" <<'PY' || fail=1
import sys
root, work = sys.argv[1:]
sys.path.insert(0, f"{root}/tests/gas")
import numpy as np
import collapse_cases as cc
import pore_cases as pc
ref = f"{work}/ref"
cc.make(ref, end=0.04)
cc.run(ref)
d = max(float(np.abs(np.asarray(pc.field(ref, "0.04", f), float)
                     - np.asarray(pc.field(f"{work}/dam", "0.04", f), float)).max())
        for f in ("alpha", "p", "U", "sigma"))
print(f"  fluidisedDamBreak against the test case at 0.04 s: largest difference {d:.1e}")
sys.exit(0 if d == 0 else 1)
PY

[ "$fail" -eq 0 ] && echo "test_gas_tutorials: PASS" || echo "test_gas_tutorials: FAIL"
exit $fail
