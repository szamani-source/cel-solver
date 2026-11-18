#!/bin/bash
# WP 4.4b: limiting lateral pressure on a pile section moving through
# undrained clay (Randolph and Houlsby 1984, smooth 9.14), the pile a moving
# STL body (decision 0053, 0056). Builds celVoidCheck and runs
# tests/pile/pile_checks.py (about 2 min).
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/pile/pile_checks.py" "$work" "$@"
