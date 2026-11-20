#!/bin/bash
# WP G.2: free surface on angled (leaning) hex cells (decision 0104): a lake
# at rest and a granular column collapse against the Cartesian mesh.
# About 2 minutes.
#
# Builds the test driver celVoidCheck and runs tests/angled/free_surface_checks.py.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT
mkdir -p "$work"

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/angled/free_surface_checks.py" "$work" "$@"
