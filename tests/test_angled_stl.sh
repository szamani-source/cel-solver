#!/bin/bash
# WP G.5: STL (cut-cell) walls with a free surface on angled cells (decision
# 0107): lakes in an STL tank and against an inclined STL wall on leaning
# cells against the Cartesian mesh. About 1 minute.
#
# Builds the test driver celVoidCheck and runs tests/angled/stl_angled_checks.py.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT
mkdir -p "$work"

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/angled/stl_angled_checks.py" "$work" "$@"
