#!/bin/bash
# WP 3.5: cohesion in the T0 tank (decision 0054). With c = 500 Pa the flow
# stops at outlets of 80 and 100 mm and runs at 120 and 140 mm (Jenike's
# critical width 109 mm); the cohesionless control flows at 80 mm. E = 1 MPa,
# 10 mm cells; the 10 MPa check is tests/manual/test_cohesion.sh.
#
# Builds the test driver tests/void/voidCheck (celVoidCheck) and runs
# tests/cohesion/cohesion_checks.py.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/cohesion/cohesion_checks.py" "$work" "$@"
