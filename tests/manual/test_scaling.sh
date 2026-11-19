#!/bin/bash
# WP 6.3a: parallel scaling of celFoam on the T0 tank (10 and 5 mm cells, 1, 2
# and 4 processors). Prints wall time, cost per cell per step and speed-up;
# results go to docs/results/scaling.md by hand. About 6 minutes on 4 cores.
set -eu
root="$(cd "$(dirname "$0")/../.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }
python3 "$root/tests/scaling/scaling_run.py" "$work" "$@"
