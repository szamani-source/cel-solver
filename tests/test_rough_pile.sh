#!/bin/bash
# WP 4.3e: rough walls (decision 0059). The pile section of WP 4.4b, fully
# rough (wall adhesion cu), against Randolph and Houlsby (1984): 11.94 within 15 % at 10 mm.
# Builds celVoidCheck and runs tests/pile/pile_checks.py --rough (about 3 min).
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/pile/pile_checks.py" "$work" --rough "$@"
