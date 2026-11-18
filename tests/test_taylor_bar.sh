#!/bin/bash
# WP 5.2: Lagrangian Taylor bar (decision 0071). A 3D quarter model of the
# copper rod of the Abaqus Benchmarks (32.4 mm, 227 m/s, 80 microseconds)
# against its C3D8R result: shortening 13.10 mm, widening 5.528 mm.
# Builds celHex8Check and runs it with -taylor (about 1.5 min).
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/hex8/hex8Check" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

celHex8Check -taylor "$@"
