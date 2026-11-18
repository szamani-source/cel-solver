#!/bin/bash
# WP 4.3c: STL cut-cell walls (decision 0047) with the free surface (void,
# decision 0031); the combination is decision 0051. Lake at rest in an STL
# tank, a block falling onto a 10 deg STL floor without leaking, and the
# WP 4.3a tilted tank with a void top against its lid version.
# Builds celVoidCheck and celSiloCheck and runs
# tests/stlwalls/stl_void_checks.py. Takes a few minutes.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

for d in void/voidCheck silo/siloCheck
do
    wmake "$root/tests/$d" > "$work/log.wmake" 2>&1 \
        || { tail -n 40 "$work/log.wmake"; exit 1; }
done

python3 "$root/tests/stlwalls/stl_void_checks.py" "$work" "$@"
