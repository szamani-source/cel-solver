#!/bin/bash
# WP T.6: a moving bed past heated walls (decision 0136), the case of Adapa
# et al. (2024). tests/thermal/moving_bed_checks.py: local Nu and mean h
# against the slug-flow solution with the near-wall gap (5 %), mean h in
# their 225-350 W/(m^2 K), energy balance 1e-12. About 1 minute.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

python3 "$root/tests/thermal/moving_bed_checks.py" "$work"
