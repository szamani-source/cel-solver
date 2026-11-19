#!/bin/bash
# Decision 0091: a rigid footing held 2.5 mm above a sand layer at rest must
# leave the layer at rest (psi 0 and psi 30). Passes since WP A.3b (decisions
# 0092, 0094: empty slivers under a body stay empty). About 1 minute.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/footing/hover_checks.py" "$work" "$@"
