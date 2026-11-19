#!/bin/bash
# WP 3.8c: cylindrical silo discharge on the wedge against 3D Beverloo and
# a 3D quarter model (decisions 0061, 0062).
#
# Builds the test driver celVoidCheck and runs
# tests/wedgesilo/wedge_silo_checks.py.
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

python3 "$root/tests/wedgesilo/wedge_silo_checks.py" "$work" "$@"
