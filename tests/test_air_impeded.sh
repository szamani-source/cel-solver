#!/bin/bash
# WP F.5: air-impeded discharge of fine powder from a plane-flow hopper:
# gas on against gas off over grain sizes 4 to 0.1 mm, 2 and 4 processors,
# cost of gas on against off (decisions 0125, 0130).
# tests/gas/air_impeded_checks.py. About 8 minutes on 4 cores.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

python3 "$root/tests/gas/air_impeded_checks.py" "$work"
