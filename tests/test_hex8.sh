#!/bin/bash
# WP 5.1: explicit hex8 Lagrangian body (decision 0070). Builds celHex8Check
# and runs it: MacNeal-Harder patch test, a 1D strain wave in a steel bar,
# energy of a free elastic cube with hourglass modes (a few seconds).
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/hex8/hex8Check" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

celHex8Check
