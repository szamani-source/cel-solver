#!/bin/bash
# A.4 done-when (decision 0095): smooth strip footing pushed into phi 30 sand
# with associated flow (psi 30), 10 mm cells, to a settlement of 0.06 B. Before
# voidPlasticFlow the partly filled surface cells beside the footing could not
# dilate, their pressure grew without bound and the run blew up. Checks: the
# run completes, the load stays bounded and steady, mass kept, no packing.
# N_gamma itself is still high (rate effect, decision 0095): manual
# tests/manual/test_sand_footing.sh. About 2 min.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/footing/sand_checks.py" "$work" --psi 30 --stable
