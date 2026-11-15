#!/bin/bash
# WP 3.3: the T0 tank with an outlet that is a closed wall while the material
# settles and opens at a set time (materialProperties outlets, decision 0033;
# prototype rules 10 and 10b). Nothing leaves before the opening, mass balance
# exact, rate at E = 0.1 MPa against the prototype (10 mm cells).
#
# Builds the test driver tests/void/voidCheck (celVoidCheck) and runs
# tests/outlet/outlet_checks.py.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/outlet/outlet_checks.py" "$work" "$@"
