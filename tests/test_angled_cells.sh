#!/bin/bash
# WP G.1: the solver on angled (non-orthogonal) hex cells (decision 0103):
# a P-wave on a 45 degree sheared mesh, and an elastic column under gravity
# on leaning cells with the solver's own time step. About 1 minute.
#
# Builds the test driver celCycleCheck and runs tests/angled/angled_checks.py.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT
mkdir -p "$work"

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/eulerian/cycleCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/angled/angled_checks.py" "$work" "$@"
