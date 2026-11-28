#!/bin/bash
# WP T.8: tutorials/solidFluidThermal/hotSiloDischarge runs (E 0.1 MPa, to 0.8 s): the energy
# balance holds to 1e-12 of the thermal energy, the sand leaving the outlet
# is the hot layer's 150 C, and no cell leaves the 20 to 150 C range of the
# start. About 30 seconds.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

cp -r "$root/tutorials/solidFluidThermal/hotSiloDischarge" "$work/"
cd "$work/hotSiloDischarge"
foamDictionary -entry endTime -set 0.8 system/controlDict > /dev/null
foamDictionary -entry writeInterval -set 0.1 system/controlDict > /dev/null
./Allrun 0.1 > /dev/null 2>&1 || { tail -20 log.celFoam; echo "FAIL: Allrun"; exit 1; }

python3 - <<'PY'
import re, sys
log = open("log.celFoam").read()
e0 = float(re.search(r"^Thermal energy at start\s*=\s*(\S+)", log, re.M).group(1))
bal = [abs(float(v)) for v in re.findall(r"^Energy balance\s*=\s*(\S+)", log, re.M)]
ok = True
rel = max(bal) / e0
print(f"energy balance {rel:.1e} of the thermal energy (limit 1e-12)")
ok &= rel < 1e-12
out = [float(v) for v in re.findall(r"^Discharged mass\s*=\s*(\S+)", log, re.M)][-1]
car = [float(v) for v in re.findall(r"^Heat carried in\s*=\s*(\S+)", log, re.M)][-1]
T_out = -car / (800.0 * out)
print(f"sand left at {T_out - 273.15:.6f} C (150 C)")
ok &= abs(T_out - 423.15) < 1e-6
s = open("0.8/T").read()
n = int(re.search(r"List<scalar>\s*(\d+)", s).group(1))
T = [float(v) for v in s[s.index("(", s.index("List<scalar>")) + 1:].split()[:n]]
print(f"T in [{min(T) - 273.15:.6f}, {max(T) - 273.15:.6f}] C")
ok &= min(T) > 293.15 - 1e-9 and max(T) < 423.15 + 1e-9
print("test_hot_silo_tutorial: " + ("PASS" if ok else "FAIL"))
sys.exit(0 if ok else 1)
PY
