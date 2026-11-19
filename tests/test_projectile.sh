#!/bin/bash
# WP 5.4: projectile into clay (target problem T3). A free slug, rigid and as
# a steel hex8 body, fired into Tresca clay: penetration against Prandtl's
# bearing resistance by energy, steel against rigid.
# Builds celVoidCheck and runs tests/projectile/projectile_checks.py (about 4 min).
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/projectile/projectile_checks.py" "$work" "$@"
