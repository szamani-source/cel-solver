#!/bin/bash
# WP 2.4: remap in C++ (eulerianSolid::remap, eulerianSolidRemap.C).
#
# Builds tests/eulerian/remapCheck (celRemapCheck), which calls remap alone
# with a prescribed motion, and runs tests/eulerian/remap_test.py:
#   A. x translation on a 64x3x3 mesh matches the 1D prototype remap
#      (lr1d.py density, lr1d_plastic.py mass-weighted fields) to 1e-10
#   B. smooth profile, one box crossing along the 2D diagonal (32/64/128) and
#      the 3D body diagonal (12/24/48): order >= 1.7 for rho, U, sigma
#   C. square profile, diagonal: no new extrema above 1e-10 in rho, U, sigma;
#      history (superbee) edge width <= 4 cells after 50 cells in x and y
#   D. slotted disc of history, one turn of solid-body rotation: total kept,
#      no new extrema above 1 %, cells without flux unchanged bit for bit
#   E. uniform state stays uniform under a non-uniform divergence-free flux,
#      and (velocity, stress, history) under a compressing flux and its reverse
#   F. case C on 2 and 4 processors matches serial to 1e-12
#   G. compression along x between walls keeps every row the same to 1e-12
#   Every case: mass, momentum and sum(m q) conserved to 1e-12 relative.
# Default material Courant number 0.25 per direction. Takes about 2 minutes.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

log="$(mktemp)"
trap 'rm -f "$log"' EXIT
wmake "$root/tests/eulerian/remapCheck" > "$log" 2>&1 \
    || { tail -n 40 "$log"; exit 1; }

python3 "$root/tests/eulerian/remap_test.py"
