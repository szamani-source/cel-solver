#!/bin/bash
# WP 3.6s: mesh dependence of the collapse shear band (decision 0038).
# Column collapse a = 0.5, E = 1 MPa, plain T0 Mohr-Coulomb and mu(I), on
# 10, 5 and 2.5 mm cells to t = 1 s (about 9 min on 4 cores). Passes when all
# six runs finish and keep mass to 1e-10; the band widths and run-outs are
# printed and recorded in docs/wp/WP-3.6s.md, not checked (they do not converge).
# Optional arguments: cell sizes in mm (default 10 5 2.5).
set -eu
root="$(cd "$(dirname "$0")/../.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/collapse/shear_band.py" "$work" "$@" | tee "$work/table"
if awk 'NR > 2 && $1 ~ /^(MC|muI)$/ { n++; if ($(NF-2) + 0 > 1e-10) bad = 1 }
        END { exit !(n >= 2 && !bad) }' "$work/table"
then
    echo "Shear-band study: PASS (all runs finished, mass kept)"
else
    echo "Shear-band study: FAIL"; exit 1
fi
