#!/bin/bash
# WP 6.2a: restart (decision 0081). celFoam runs the T0 tank (outlet open,
# discharging) and a tank with a free block and a hex8 footing in contact,
# continuous and restarted from a written time (serial, and the bodies case on
# 2 processors): fields, body state, mass counters and bodyForces.dat at the
# end identical bit for bit. Needs celFoam built (about 40 s).
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }
command -v celFoam > /dev/null || { echo "FAIL: celFoam not built (./Allwmake)"; exit 1; }

python3 "$root/tests/restart/restart_checks.py" "$work" "$@"
