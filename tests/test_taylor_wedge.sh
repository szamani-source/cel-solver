#!/bin/bash
# WP 2.6: the Taylor bar (copper rod at 227 m/s) on the axisymmetric wedge
# against the axisymmetric Abaqus benchmark (decisions 0061, 0063).
#
# Builds the test driver celVoidCheck and runs
# tests/taylorwedge/taylor_wedge_checks.py.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT
mkdir -p "$work"

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

for d in void/voidCheck; do
    wmake "$root/tests/$d" > "$work/log.wmake" 2>&1 \
        || { tail -n 40 "$work/log.wmake"; exit 1; }
done

python3 "$root/tests/taylorwedge/taylor_wedge_checks.py" "$work" "$@"
