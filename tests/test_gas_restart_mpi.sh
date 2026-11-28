#!/bin/bash
# WP F.9: the gas runs restart and run in parallel (decision 0141).
# tests/gas/restart_mpi_checks.py: the F.5 hopper and the F.8 dam break
# restarted match the continuous runs to 1e-12 (every field); the hopper on
# 2 and 4 processors matches serial to 1e-10; the dam break with cuts through
# the column is reported (it parts when round-off flips a free-surface rule).
# About 4 minutes.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

python3 "$root/tests/gas/restart_mpi_checks.py" "$work"
