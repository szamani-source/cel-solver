#!/bin/bash
# WP 3.4: T0 tank discharge in celFoam with the collapse code (shear cap,
# rule 9) at E = 0.1 MPa: mass balance, no stress in void, and two processors
# (outlet split across them) against serial, rate within 2 percent (decision
# 0034). The E sweep is tests/manual/test_discharge_stiffness.sh.
#
# Builds the test driver tests/void/voidCheck (celVoidCheck) and runs
# tests/discharge/discharge_checks.py.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/discharge/discharge_checks.py" "$work" "$@"
