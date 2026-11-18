#!/bin/bash
# WP 4.5: silo and hopper discharge with STL walls (decisions 0047, 0051,
# 0052). The T0 tank with STL walls and floor against mesh walls (bands of
# three runs), and a 60 deg STL hopper run to empty (rate recorded).
# Builds celVoidCheck and runs tests/stlwalls/stl_silo_checks.py. Takes
# about 5 minutes.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/stlwalls/stl_silo_checks.py" "$work" "$@"
