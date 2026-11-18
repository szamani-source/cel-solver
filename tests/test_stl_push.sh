#!/bin/bash
# WP 4.3b: moving STL bodies (decision 0053). A footing pushed into a T0
# layer and pulled out again: no mass lost, no gap under it, momentum
# balance with the reported body force, separation on pull-out.
# Builds celVoidCheck and runs tests/stlwalls/stl_push_checks.py (about 1 min).
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/stlwalls/stl_push_checks.py" "$work" "$@"
