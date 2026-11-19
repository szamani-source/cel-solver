#!/bin/bash
# WP 3.7b: the T0 tank discharging while an inlet in its lid feeds it
# (inlets, decision 0060). Discharge rate independent of the
# feed, a matched feed holds the level, exact mass balance.
#
# Builds the test driver tests/void/voidCheck (celVoidCheck) and runs
# tests/fedsilo/fed_silo_checks.py.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/fedsilo/fed_silo_checks.py" "$work" "$@"
