#!/bin/bash
# WP 4.4: strip footing on undrained clay (Prandtl, Nc = 5.14), the footing a
# moving STL body (decision 0053). Builds celVoidCheck and runs
# tests/footing/footing_checks.py (about 1 min).
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/footing/footing_checks.py" "$work" "$@"
