#!/bin/bash
# WP K.4: compaction cap with STL (immersed) walls (decision 0122): the K.2
# settling column held between frictionless STL walls lies on the line
# and settles as the mesh-wall column within 1 %.
# About 1 minute.
#
# Builds the test driver celVoidCheck and runs tests/compaction/settle_column_stl.py.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT
mkdir -p "$work"

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/compaction/settle_column_stl.py" "$work" "$@"
