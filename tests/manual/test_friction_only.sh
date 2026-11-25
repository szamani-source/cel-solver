#!/bin/bash
# WP F.8g: pore pressure that only lowers the friction (gasProperties
# frictionOnly; decision 0133). tests/gas/friction_only_checks.py:
# full support spreads as frictionless, support lost, 2 processors.
set -eu
root="$(cd "$(dirname "$0")/../.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

python3 "$root/tests/gas/friction_only_checks.py" "$work"
