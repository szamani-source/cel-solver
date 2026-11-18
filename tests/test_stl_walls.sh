#!/bin/bash
# WP 4.3a: fixed STL walls with friction as cut-cell walls (decisions 0042,
# 0047). The WP 3.1 frictional tank with its walls given as STL boxes,
# aligned with the mesh and turned 10 deg against it, against the mesh-wall
# tank: Janssen stress profile, leak, force balance.
# Builds tests/silo/siloCheck (celSiloCheck) and runs
# tests/stlwalls/stl_wall_checks.py. Takes a few minutes.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/silo/siloCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/stlwalls/stl_wall_checks.py" "$work"
