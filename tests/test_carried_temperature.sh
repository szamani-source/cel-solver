#!/bin/bash
# WP T.2: temperature moves with the material (decision 0135).
# tests/thermal/carried_checks.py: a hot-cold block carried 50 cells through
# void (energy 1e-12, front sharp and in place), a box fed through an inlet and
# emptied through an outlet (energy balance 1e-12, fed material at the inlet
# temperature, 2 processors against serial 1e-10), one temperature stays exact.
# About 2-3 minutes.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

python3 "$root/tests/thermal/carried_checks.py" "$work"
