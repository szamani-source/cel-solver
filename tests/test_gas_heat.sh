#!/bin/bash
# WP T.4: gas temperature and grain-gas exchange (decision 0136).
# tests/thermal/gas_heat_checks.py: hot air into a cold packed column against
# Schumann (1929) with Gunn's (1978) h (1 %), equilibrium front speed (1 %),
# uniform temperature under flow (1e-12), energy balance, 2 processors.
# About 1 minute.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

python3 "$root/tests/thermal/gas_heat_checks.py" "$work"
