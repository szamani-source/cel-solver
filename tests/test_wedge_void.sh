#!/bin/bash
# WP 3.8b: the free surface on an axisymmetric wedge mesh (decisions 0031,
# 0061): lake at rest, falling plug, void-top tank, 2 processors.
#
# Builds the test drivers celVoidCheck and celSiloCheck and runs
# tests/wedgevoid/wedge_void_checks.py.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT
mkdir -p "$work"

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

for d in void/voidCheck silo/siloCheck; do
    wmake "$root/tests/$d" > "$work/log.wmake" 2>&1 \
        || { tail -n 40 "$work/log.wmake"; exit 1; }
done

python3 "$root/tests/wedgevoid/wedge_void_checks.py" "$work" "$@"
