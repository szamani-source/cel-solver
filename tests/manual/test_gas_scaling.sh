#!/bin/bash
# WP F.9: cost of the gas and parallel speed-up (decision 0141).
# tests/gas/scaling_checks.py: F.5 hopper refined to 8,704 cells, 200 steps;
# gas on at most 1.5 times gas off per step, 4 processors at least 2.5 times
# serial. Needs 4 idle cores. About 3 minutes.
set -eu
root="$(cd "$(dirname "$0")/../.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

python3 "$root/tests/gas/scaling_checks.py" "$work" "${1:-200}"
