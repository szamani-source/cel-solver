#!/bin/bash
# WP F.9a (decision 0134): after the fluidised bed lifts, a serial run with
# only pTolerance changed parts from serial as fast as 2 and 4 processors
# do: the bed amplifies round-off, the parallel code is not at fault.
# tests/gas/parallel_checks.py --lift. Report only. About 4 minutes.
set -eu
root="$(cd "$(dirname "$0")/../.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

python3 "$root/tests/gas/parallel_checks.py" "$work" --lift
