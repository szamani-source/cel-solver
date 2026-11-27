#!/bin/bash
# WP T.0: temperature reference values. Runs the self-checks of
# tools/thermal_reference.py (air properties against Incropera Table A.4, erfc and the
# convectively cooled slab against finite differences and Incropera Table 5.1,
# Zehner-Schluender limits and bounds, Gunn limits, Schumann against an independent
# finite-difference solution and the equilibrium front speed, Graetz plates against
# Shah and London, slug flow with a wall gap against marching) and checks that the
# generated block in docs/spec/thermal-references.md is exactly what the script
# prints. Pure Python, no OpenFOAM needed (about 20 s).
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
tool="$root/tools/thermal_reference.py"
doc="$root/docs/spec/thermal-references.md"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

python3 "$tool" --check

python3 "$tool" > "$tmp/expected"
awk '/<!-- BEGIN thermal_reference.py output -->/{on=1} on{print} /<!-- END thermal_reference.py output -->/{on=0}' \
    "$doc" > "$tmp/actual"
if diff -u "$tmp/expected" "$tmp/actual"; then
    echo "ok  table in docs/spec/thermal-references.md matches tools/thermal_reference.py"
else
    echo "BAD table in docs/spec/thermal-references.md differs from tools/thermal_reference.py"
    exit 1
fi
echo "test_thermal_reference: PASS"
