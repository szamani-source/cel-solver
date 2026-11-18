#!/bin/bash
# WP 3.6g: dilatancy in the mu(I) Mohr-Coulomb and Drucker-Prager returns
# (decision 0048). Kernel checks in tests/constitutive/dilatancyCheck.C, a
# plane-strain shear box at constant normal load (about 2 s):
# dense start dilates and loose start contracts to phiEq(I) = phiC - bPhi I
# (1e-6), steady friction equals plain mu(I), flow rule sin(psi) at the first
# plastic step, off without the entries, Drucker-Prager reaches phiEq too.
# Bit-for-bit equality of the models without the entries is checked by
# tests/test_mu_i.sh and tests/test_constitutive_cpp.sh.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
kernel="$root/src/celConstitutive/kernel"

g++ -std=c++17 -O2 -Wall -Wextra -Werror -I"$kernel" \
    "$root/tests/constitutive/dilatancyCheck.C" "$kernel"/*.C \
    -o "$work/dilatancyCheck"
"$work/dilatancyCheck"
echo PASS
