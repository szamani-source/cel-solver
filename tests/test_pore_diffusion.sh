#!/bin/bash
# WP F.8e: pore pressure as an advected-diffused field (gasProperties
# diffusivity; decision 0133). tests/gas/pore_diffusion_checks.py: step at a
# face against erfc with the given diffusivity at two values, 2 and 4
# processors against serial.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

python3 "$root/tests/gas/pore_diffusion_checks.py" "$work"
