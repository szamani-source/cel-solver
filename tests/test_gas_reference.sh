#!/bin/bash
# WP F.0: gas reference values for version 2. Runs the self-checks of
# tools/gas_reference.py (Ergun -> Darcy, Terzaghi and erfc against finite
# differences, Wen-Yu constants from Ergun, Di Felice and Richardson-Zaki limits,
# Stokes settling) and checks that the generated block in docs/spec/gas-references.md
# is exactly what the script prints. Pure Python, no OpenFOAM needed.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
tool="$root/tools/gas_reference.py"
doc="$root/docs/spec/gas-references.md"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

python3 "$tool" --check

python3 "$tool" > "$tmp/expected"
awk '/<!-- BEGIN gas_reference.py output -->/{on=1} on{print} /<!-- END gas_reference.py output -->/{on=0}' \
    "$doc" > "$tmp/actual"
if diff -u "$tmp/expected" "$tmp/actual"; then
    echo "ok  table in docs/spec/gas-references.md matches tools/gas_reference.py"
else
    echo "BAD table in docs/spec/gas-references.md differs from tools/gas_reference.py"
    exit 1
fi
echo "test_gas_reference: PASS"
