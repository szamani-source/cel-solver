#!/bin/bash
# WP F.6: gas momentum in the void (model gasFlow): plane channel against
# Poiseuille, still air in a closed box, time step unchanged, bed and void
# together (venting column, outlet box), 2 and 4 processors, cost of gas
# on against off (decisions 0125, 0131). tests/gas/gas_flow_checks.py.
# About 3 minutes.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

python3 "$root/tests/gas/gas_flow_checks.py" "$work"
