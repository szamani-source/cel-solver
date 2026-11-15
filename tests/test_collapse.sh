#!/bin/bash
# WP 3.2b: granular column collapse in celFoam's cycle with the face shear cap
# (rule 8) and rule 9, neither rule 12 nor 13 (decisions 0023, 0024, 0029,
# 0032). a = 0.5 at E = 0.1, 1, 10 MPa and a = 2 at 1 MPa, 10 mm cells; the
# 1 MPa runs are checked against decision 0019, the others recorded.
#
# Builds the test driver tests/void/voidCheck (celVoidCheck) and runs
# tests/collapse/collapse_checks.py.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/collapse/collapse_checks.py" "$work" "$@"
