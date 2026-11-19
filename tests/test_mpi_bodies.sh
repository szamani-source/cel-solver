#!/bin/bash
# WP 6.1: MPI runs with Lagrangian bodies (decision 0080). A free block
# dropped onto T0 and a hex8 footing pushed into it (2 processors split
# through the body), and the T0 tank with STL walls and an outlet (4
# processors, 2 x 2): first 100 steps match serial to 1e-12, outcomes within
# the round-off band of serial (a perturbed serial run is printed beside).
# Builds celVoidCheck and runs tests/mpi/mpi_checks.py (about 8 min).
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/mpi/mpi_checks.py" "$work" "$@"
