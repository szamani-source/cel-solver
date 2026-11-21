#!/bin/bash
# Top-fill fix step 2 (decision 0114): rest relaxation removes the
# impact-locked lateral stress of a box filled from the top (K 0.87 ->
# about K0) and leaves a box placed with gravity on unchanged.
# About 6 minutes.
#
# Builds the test driver celVoidCheck and runs tests/compaction/relax_checks.py.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT
mkdir -p "$work"

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/compaction/relax_checks.py" "$work" "$@"
