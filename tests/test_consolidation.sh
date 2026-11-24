#!/bin/bash
# WP F.3: gas-filled column loaded at the top against Terzaghi, column at
# rest, 2 and 4 processors, cost of gas on against off (decisions 0125,
# 0128). tests/gas/consolidation_checks.py. About 2 minutes.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

python3 "$root/tests/gas/consolidation_checks.py" "$work"
