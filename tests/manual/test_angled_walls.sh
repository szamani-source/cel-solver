#!/bin/bash
# WP G.3: walls and free surface on angled faces (decision 0106): an elastic
# layer sliding down, or resting on, a 30 degree frictional slope on a
# Cartesian, a turned and a sheared mesh. About 2 minutes.
# In tests/manual until G.3 passes (decision 0106); then it moves to tests/.
#
# Builds the test driver celVoidCheck and runs tests/angled/slope_checks.py.
set -eu
root="$(cd "$(dirname "$0")/../.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT
mkdir -p "$work"

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/angled/slope_checks.py" "$work" "$@"
