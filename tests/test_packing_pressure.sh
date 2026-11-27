#!/bin/bash
# WP F.8f: packing-dependent frictional pressure (materialProperties
# princetonFr; decision 0133). tests/gas/packing_checks.py: the bound
# P <= Pc(es), compaction of a column, option off, 2 processors.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

python3 "$root/tests/gas/packing_checks.py" "$work"
