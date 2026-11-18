#!/bin/bash
# WP 3.7a: a solid inlet (materialProperties inlets, decision 0056) feeds
# fresh material into an empty box: closed before its start time, exact mass
# rate while open, the material lands and comes to rest at the inlet density.
#
# Builds the test driver tests/void/voidCheck (celVoidCheck) and runs
# tests/inlet/inlet_checks.py.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/inlet/inlet_checks.py" "$work" "$@"
