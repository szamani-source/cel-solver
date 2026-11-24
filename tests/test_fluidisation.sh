#!/bin/bash
# WP F.7: gas through and around the grains: pore gas momentum against its
# closed form, onset of fluidisation against Wen and Yu, dilute settling
# against Richardson and Zaki, Stage A reproduced, 2 and 4 processors, cost
# (decisions 0125, 0132). tests/gas/fluidisation_checks.py. About 5 minutes.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

python3 "$root/tests/gas/fluidisation_checks.py" "$work"
