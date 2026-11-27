#!/bin/bash
# WP T.7 done-when (decision 0136): gas density, viscosity and conductivity
# from the gas temperature (thermalProperties variableGas) and plastic
# heating (plasticHeating). Sealed box p/T to 1e-10, gas mass 1e-12,
# 2 processors 1e-10, plastic heating against the closed form 1e-10, energy
# balance with the plastic heat 1e-12, T.4 and T.5 within their limits with
# variableGas. About 5 minutes. See tests/thermal/gas_state_checks.py.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
python3 "$root/tests/thermal/gas_state_checks.py" "$work"
