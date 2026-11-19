#!/bin/bash
# WP 4.4c: smooth strip footing on sand (phi 30) against the exact N_gamma 7.653
# of Martin (2005). Not passing yet (decision 0057): with psi 0 the footing reads
# N_gamma about 11 at 10 mm and 9.7 at 5 mm; with psi 30 (default) the run blows up.
# Options pass through: --fine (5 mm cells, about 11 min), --psi DEG, and the
# options of sand_checks.py. Decision 0095: the reading depends on the push speed
# times rho c (numerical viscosity); psi 30 now runs (voidPlasticFlow) and reads
# about 20 at the default speed, 11.9 at an eighth of it.
set -eu
root="$(cd "$(dirname "$0")/../.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/footing/sand_checks.py" "$work" "$@"
