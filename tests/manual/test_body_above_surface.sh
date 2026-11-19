#!/bin/bash
# Decision 0091: a rigid footing held 2.5 mm above a sand layer at rest must
# leave the layer at rest (psi 0 and psi 30). Not passing yet: it is the
# done-when of the next package (body near a free surface). About 1 minute.
set -eu
root="$(cd "$(dirname "$0")/../.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/footing/hover_checks.py" "$work" "$@"
