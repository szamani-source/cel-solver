#!/bin/bash
# WP H.2: the conical hopper wall on the axisymmetric wedge (decision 0101):
# the Couto silo at rest with a cut-cell hopper wall. About 4 minutes.
#
# Builds the test driver celVoidCheck and runs tests/hopper/hopper_wedge_checks.py.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT
mkdir -p "$work"

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/hopper/hopper_wedge_checks.py" "$work" "$@"
