#!/bin/bash
# WP 3.0: silo reference values. Runs the self-checks of tools/silo_reference.py
# (wall-yield K closed form vs Mohr circle, Janssen ODE vs closed form, force balance,
# Beverloo scaling) and checks that the tables in docs/spec/silo-reference.md are
# exactly what the script prints. Pure Python, no OpenFOAM needed.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
tool="$root/tools/silo_reference.py"
doc="$root/docs/spec/silo-reference.md"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

python3 "$tool" --check

# Expected: the two generated blocks, in the order they appear in the page.
{ python3 "$tool"; python3 "$tool" --H 1.2; } > "$tmp/expected"
awk '/<!-- BEGIN silo_reference.py output -->/{on=1} on{print} /<!-- END silo_reference.py output -->/{on=0}' \
    "$doc" > "$tmp/actual"

if diff -u "$tmp/expected" "$tmp/actual"; then
    echo "ok  tables in docs/spec/silo-reference.md match tools/silo_reference.py"
else
    echo "BAD tables in docs/spec/silo-reference.md differ from tools/silo_reference.py"
    exit 1
fi

# Spot values a reader would quote (guards against a silently changed formula).
python3 - "$tool" <<'PY'
import importlib.util, math, sys
spec = importlib.util.spec_from_file_location("sr", sys.argv[1])
sr = importlib.util.module_from_spec(spec); spec.loader.exec_module(sr)
phi, dw = math.radians(35.0), math.radians(20.0)
checks = [
    ("Rankine active phi 35", sr.k_rankine_active(phi), 0.270990, 1e-6),
    ("Jaky phi 35", sr.k_jaky(phi), 0.426424, 1e-6),
    ("elastic nu 0.3", sr.k_elastic(0.3), 0.428571, 1e-6),
    ("wall yield phi 35 delta 20", sr.k_wall_yield(phi, dw), 0.290608, 1e-6),
    ("Beverloo D 0.06 C 1.45", sr.beverloo(1.45, 1800.0, 9.81, 0.06), 120.14, 0.01),
]
bad = 0
for name, got, want, tol in checks:
    ok = abs(got - want) <= tol
    bad += not ok
    print(f"{'ok ' if ok else 'BAD'} {name}: {got:.6f} (expected {want})")
sys.exit(1 if bad else 0)
PY
