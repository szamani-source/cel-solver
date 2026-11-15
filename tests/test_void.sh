#!/bin/bash
# WP 3.2a: one material plus void in celFoam's cycle (free surface, decision 0031),
# checked on the prototype's cases (WP 3.2p, decision 0021): lake at rest with a
# half-full top cell, block translating through void, free fall and landing,
# cohesionless block pulled apart, 2 processors against serial.
#
# Builds the test driver tests/void/voidCheck (celVoidCheck) and runs
# tests/void/void_checks.py.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/void/void_checks.py" "$work" "$@"
