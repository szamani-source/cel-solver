#!/bin/bash
# WP 5.3a: two-way coupling of deformable hex8 bodies (decision 0072). The
# footing push of WP 4.3b as a rigid and a steel body: the steel footing
# follows the rigid force, mass and density kept.
# Builds celVoidCheck and runs tests/deformable/deformable_push_checks.py (about 3.5 min).
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/deformable/deformable_push_checks.py" "$work" "$@"
