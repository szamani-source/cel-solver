#!/bin/bash
# WP T.3: convective and heat-flux walls, Zehner-Schluender bed conductivity
# (decision 0136). tests/thermal/walls_checks.py: slab cooled at Bi 1 and 10
# against the Carslaw-Jaeger series (1 %), energy balance 1e-12, 2 processors
# 1e-10, heat flux q A t exact, kBed against tools/thermal_reference.py 1e-12.
# About 1 minute.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

python3 "$root/tests/thermal/walls_checks.py" "$work"
