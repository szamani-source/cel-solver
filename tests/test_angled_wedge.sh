#!/bin/bash
# WP G.4: angled cells on the axisymmetric wedge: the Couto silo at rest on
# the wall-aligned mesh (free surface, T0 rules); wall pressures within 5 %
# of Janssen with K = nu/(1-nu). About 13 minutes.
#
# Builds the test driver celVoidCheck and runs tests/angled/wedge_aligned_checks.py.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT
mkdir -p "$work"

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/angled/wedge_aligned_checks.py" "$work" "$@"
