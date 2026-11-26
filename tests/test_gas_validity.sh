#!/bin/bash
# Range of validity of the gas coupling (decision 0135). The lift check
# (gasProperties liftCheck, on by default) warns when the gas carries the
# weight of the grains: on the fluidisation bed it warns at the onset
# (within 10 % of U_mf), not before; on the fine-powder hopper with the
# outlet opening (inside the range, validation pack test 4) it does not
# warn. tests/gas/validity_checks.py. About 2 minutes.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

python3 "$root/tests/gas/validity_checks.py" "$work"
