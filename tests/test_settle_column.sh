#!/bin/bash
# WP K.2: compaction cap in the Eulerian solver (decision 0113): a column
# settling under its own weight lies on the compression line
# ln(rho/rho0) = P0/K + lambda* ln(P/P0) cell by cell, mass conserved.
# About 1 minute.
#
# Builds the test driver celVoidCheck and runs tests/compaction/settle_column.py.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT
mkdir -p "$work"

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/compaction/settle_column.py" "$work" "$@"
