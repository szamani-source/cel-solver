#!/bin/bash
# Regression tests. Run from anywhere with an OpenFOAM v2512 environment loaded.
# Exit code is non-zero if any test fails.
set -u
root="$(cd "$(dirname "$0")/.." && pwd)"
fail=0

for t in "$root"/tests/test_*.sh; do
    name="$(basename "$t" .sh)"
    if bash "$t" > "$root/tests/$name.log" 2>&1; then
        echo "PASS  $name"
    else
        echo "FAIL  $name  (see tests/$name.log)"
        fail=1
    fi
done

exit $fail
