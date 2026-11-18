#!/bin/bash
# WP 3.5 (decision 0054), beyond CI: the cohesion bracket at T0's E = 10 MPa
# (c = 500 Pa: 100 mm stops, 120 mm flows) and the scaling with c (c = 250 Pa:
# Jenike 55 mm; 40 mm creeps below 5 % of Beverloo, 80 mm flows). Prints one
# line per run; judge with the stop / flow thresholds of
# tests/cohesion/cohesion_checks.py. About 2 minutes on 4 cores.
set -eu
root="$(cd "$(dirname "$0")/../.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

chk="$root/tests/cohesion/cohesion_checks.py"
for args in "500 0.10 10" "500 0.12 10" "250 0.04 1" "250 0.08 1"
do
    # shellcheck disable=SC2086
    python3 "$chk" "$work" --case $args > "$work/out.${args// /_}" 2>&1 &
done
wait
cat "$work"/out.*
