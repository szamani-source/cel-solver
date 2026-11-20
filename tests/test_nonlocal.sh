#!/bin/bash
# WP N.1: nonlocal granular fluidity (decision 0111): the kernel's point
# update gives the local law mu_s + b I in homogeneous shear, and the
# fluidity field across a layer with mu through mu_s matches an independent
# fine-grid solution (full layer and with void on top). Under a minute.
#
# Builds the test driver celFluidityCheck and runs tests/nonlocal/nonlocal_checks.py.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT
mkdir -p "$work"

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/nonlocal/fluidityCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/nonlocal/nonlocal_checks.py" "$work" "$@"
