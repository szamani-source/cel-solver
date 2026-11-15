#!/bin/bash
# WP 3.1: celFoam's cycle settles a tank under gravity with smooth and frictional
# walls and matches the prototype (WP 3.1p, decisions 0020 and 0030).
#
# Builds the test driver tests/silo/siloCheck (celSiloCheck) and runs
# tests/silo/silo_checks.py: lake at rest, smooth-wall column (K0, rho g H),
# frictional tank against Janssen (wall ratio, force balance), 2 processors
# against serial. Takes about a minute.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/silo/siloCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/silo/silo_checks.py" "$work"
