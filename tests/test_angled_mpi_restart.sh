#!/bin/bash
# WP G.6: MPI and restart on angled (leaning) cells (decision 0110): a column
# collapse on 4 processors against serial, a discharging leaning-block silo
# restarted against the continuous run (serial and 2 processors). About 1 minute.
#
# Builds the test driver celVoidCheck and runs tests/angled/mpi_restart_checks.py.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT
mkdir -p "$work"

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

command -v celFoam > /dev/null || { echo "FAIL: celFoam not built (./Allwmake)"; exit 1; }

python3 "$root/tests/angled/mpi_restart_checks.py" "$work" "$@"
