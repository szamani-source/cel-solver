#!/bin/bash
# WP T.5: gas heat in the void (decision 0135). tests/thermal/void_heat_checks.py:
# laminar air between plates at a fixed wall temperature (model gasFlow,
# twoTemperature) against the developed Nu 7.5407 (3 %), steady profile,
# energy balance 1e-10. About 2 minutes.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

python3 "$root/tests/thermal/void_heat_checks.py" "$work"
