#!/bin/bash
# WP 3.8a: axisymmetric runs on a wedge mesh (decision 0061): Lame thick
# cylinder, lake at rest and Janssen tank in a cylinder, 2 processors.
#
# Builds the test driver tests/silo/siloCheck (celSiloCheck) and runs
# tests/wedge/wedge_checks.py.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT
mkdir -p "$work"

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/silo/siloCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/wedge/wedge_checks.py" "$work" "$@"
