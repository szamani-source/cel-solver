#!/bin/bash
# WP 3.6e done-when (decision 0050): granular column collapse with the T0
# rule set mu(I) + flowImpedance + normalCap (decision 0037 placeholders),
# 10 mm cells, a = 0.5 and 2 at E = 1 MPa and T0's E = 10 MPa. Every case is
# checked against decision 0019: run-out 0.30 to 0.90 (a = 0.5) and 1.2 to
# 3.6 (a = 2), top kept within a cell at a = 0.5, fitted slope <= 37 deg,
# mass to 1e-10, no stress in empty cells. Four cases in parallel, about
# 2.5 minutes on 4 cores.
#
#   tests/manual/test_collapse_mui.sh                # all four cases
#   tests/manual/test_collapse_mui.sh 0.5:10         # a:E_MPa
set -eu
root="$(cd "$(dirname "$0")/../.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

export CEL_EXTRA_MATERIAL="mu2 0.961813;
I0 0.279;
d 0.001;
rhoS 2650.0;
flowImpedance yes;
normalCap yes;
${CEL_EXTRA_MATERIAL:-}"
export CEL_COLLAPSE_CHECK_ALL=1

cases=("$@")
[ ${#cases[@]} -gt 0 ] || cases=(0.5:10 2:10 0.5:1 2:1)
pids=()
for c in "${cases[@]}"
do
    mkdir -p "$work/$c"
    python3 "$root/tests/collapse/collapse_checks.py" "$work/$c" "$c" \
        > "$work/$c.log" 2>&1 &
    pids+=($!)
done
fail=0
for i in "${!pids[@]}"
do
    wait "${pids[$i]}" || fail=1
    grep -E "a = |run-out|checked|recorded" "$work/${cases[$i]}.log"
done
if [ $fail -eq 0 ]; then echo "Collapse with mu(I): ALL PASS"; else echo "Collapse with mu(I): FAILED"; fi
exit $fail
