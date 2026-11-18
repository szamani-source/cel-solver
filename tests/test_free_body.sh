#!/bin/bash
# WP 4.3d: free STL bodies (decision 0058). A block thrown in the void
# follows its parabola; a steel block dropped onto a T0 layer comes to rest
# with its weight carried, momentum of block and material kept.
# Builds celVoidCheck and runs tests/freebody/free_body_checks.py (about 2 min).
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/freebody/free_body_checks.py" "$work" "$@"
