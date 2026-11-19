#!/bin/bash
# WP H.3 (manual): Couto silo wall pressure at rest and in discharge against the
# experiment (decisions 0100 to 0102). About 12 minutes per model at 25 mm cells.
#   tests/manual/test_hopper_pressure.sh [mc|dp]
#   CEL_CELL=0.0125 CEL_NPROCS=4 tests/manual/test_hopper_pressure.sh mc
set -eu
root="$(cd "$(dirname "$0")/../.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT
mkdir -p "$work"
[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }
wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }
python3 "$root/tests/hopper/hopper_pressure_checks.py" "$work" "$@"
