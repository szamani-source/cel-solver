#!/bin/bash
# WP T.8: temperature restarts and runs in parallel (decision 0136).
# tests/thermal/thermal_restart_mpi_checks.py: the T.4 Schumann column and
# the T.6 moving bed restarted and on 2 and 4 processors match serial to
# 1e-10; temperature costs at most 1.2 times a gas-on step. About 3 minutes.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

python3 "$root/tests/thermal/thermal_restart_mpi_checks.py" "$work"
