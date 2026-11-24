#!/bin/bash
# WP F.4: gas with a free surface, void and an outlet: venting column
# against Terzaghi, outlet sealed until it opens and at p0 after, 2 and 4
# processors, cost of gas on against off (decisions 0125, 0129).
# tests/gas/venting_checks.py. About 3 minutes.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

python3 "$root/tests/gas/venting_checks.py" "$work"
