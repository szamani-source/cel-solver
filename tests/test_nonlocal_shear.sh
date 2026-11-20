#!/bin/bash
# WP N.2: a granular layer sheared by a moving lid under gravity, with the
# nonlocal fluidity law (decision 0111): the band width converges on 1 and
# 0.5 d cells and the velocity profile matches the 1D nonlocal solution
# within 5 % of the lid speed. About 3 minutes.
#
# Runs tests/nonlocal/shear_checks.py.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT
mkdir -p "$work"

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

python3 "$root/tests/nonlocal/shear_checks.py" "$work" "$@"
