#!/bin/bash
# WP F.7: the gas-grain drag laws ergunWenYu, ergun, wenYu, diFelice,
# syamlalOBrien and beetstra against their closed forms in tools/gas_reference.py
# (1e-10), their limits, and a permeameter with each (decisions 0125, 0132).
# tests/gas/drag_law_checks.py. About 1 minute.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

python3 "$root/tests/gas/drag_law_checks.py" "$work"
