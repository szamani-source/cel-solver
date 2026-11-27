#!/bin/bash
# WP F.9a: parallel runs with gas (decision 0134). Processor boundaries of
# every gas and solid field consistent to round-off (gasProperties
# parallelCheck) on a fluidised bed (2 x 2, 1 x 4) and a hopper with gas
# flow (2 x 2); serial and parallel agree before lift to the round-off of
# the absolute pressure. tests/gas/parallel_checks.py. About 4 minutes.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

python3 "$root/tests/gas/parallel_checks.py" "$work"
