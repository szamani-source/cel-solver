#!/bin/bash
# WP K.1: compaction cap at one point (decision 0113): isotropic and
# oedometer loading follow eps_v = lambda* ln(P/P0), unloading is elastic,
# Mohr-Coulomb and Drucker-Prager. Seconds.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT
mkdir -p "$work"

python3 "$root/tests/compaction/cap_checks.py" "$work" "$@"
