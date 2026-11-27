"""WP T.2b done-when: temperature with immersed STL walls (decision 0136).

tutorials/stlHopper (the T0 tank with two STL hopper walls at 60 deg down to
the 0.06 m outlet; E 0.1 MPa) with thermalProperties model equilibrium,
c 800 J/(kg K), k 0.25 W/(m K); the STL walls are adiabatic. The outlet opens
at 0.3 s; t = 1.5 s (the hopper is a third empty).

A  A hot lower layer (y < 0.3 m at 150 C) under 20 C: energy balance (thermal
   energy - start - heat in - heat carried in) under 1e-12 of the thermal
   energy at every write; every cell with material between 20 and 150 C (the
   small cut cells are mixed by mass and the covered cells hold no heat).
B  One temperature everywhere: T stays 20 C in every cell with material to
   1e-12, and the heat carried out is c T times the discharged mass to 1e-12.

    python3 stl_heat_checks.py WORKDIR
"""

import os
import re
import shutil
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
import conduction_checks as cc  # noqa: E402
import carried_checks as ca  # noqa: E402

report = cc.report
failed = cc.failed
T_COLD, T_HOT, C_SAND = 293.15, 423.15, 800.0
END = 1.5


def make(path, layered):
    shutil.copytree(os.path.join(ROOT, "tutorials", "stlHopper"), path)
    cd = f"{path}/system/controlDict"
    txt = open(cd).read()
    txt = re.sub(r"^endTime\s+[^;]+;", f"endTime         {END!r};", txt, flags=re.M)
    txt = re.sub(r"^writeInterval\s+[^;]+;", "writeInterval   0.25;", txt, flags=re.M)
    txt = re.sub(r"^writePrecision\s+[^;]+;", "writePrecision  17;", txt, flags=re.M)
    open(cd, "w").write(txt)
    cc.write(f"{path}/0.orig/T", cc.HEADER.format(cls="volScalarField", obj="T") + f"""
dimensions [0 0 0 1 0 0 0];
internalField uniform {T_COLD!r};
boundaryField {{ ".*" {{ type zeroGradient; }} frontAndBack {{ type empty; }} }}
""")
    if layered:
        sf = f"{path}/system/setFieldsDict"
        txt = open(sf).read()
        i = txt.rindex(");")
        txt = txt[:i] + f"""    boxToCell
    {{
        box (-1 -1 -1) (1 0.3 1);
        fieldValues ( volScalarFieldValue T {T_HOT!r} );
    }}
""" + txt[i:]
        open(sf, "w").write(txt)
    cc.write(f"{path}/constant/thermalProperties",
             cc.HEADER.format(cls="dictionary", obj="thermalProperties") + f"""
model        equilibrium;
c            {C_SAND!r};
conductivity constant;
k            0.25;
voidT        {T_COLD!r};
""")
    cc.sh(path, ["./Allrun"], "Allrun")


def material_T(path, t):
    T = cc.read_internal(f"{path}/{t:g}/T")
    rho = cc.read_internal(f"{path}/{t:g}/rho")
    if T is None:
        T = np.full(rho.size, float(open(f"{path}/{t:g}/T").read().split("internalField")[1].split()[1].rstrip(";")))
    return T[rho > 0]


def check(work):
    print("[A] Hot layer under cold sand in an STL hopper", flush=True)
    p = os.path.join(work, "layered")
    t0 = time.perf_counter()
    make(p, True)
    secs = time.perf_counter() - t0
    log = ca.energy_log(p)
    rel = max(abs(b) for b in log["bal"]) / log["E"][0]
    report("energy balance with STL walls and an outlet", rel < 1e-12,
           f"largest |balance| {rel:.1e} of the thermal energy; discharged {log['dis'][-1]:.3f} kg, "
           f"carried in {log['carried'][-1]:.1f} J, {secs:.0f} s")
    lo, hi = 1e30, -1e30
    for t in (0.5, 1.0, 1.5):
        Tm = material_T(p, t)
        lo, hi = min(lo, Tm.min()), max(hi, Tm.max())
    report("material stays between 20 and 150 C", lo > T_COLD - 1e-9 and hi < T_HOT + 1e-9,
           f"T in [{lo - 273.15:.9f}, {hi - 273.15:.9f}] C")

    print("[B] One temperature everywhere", flush=True)
    p = os.path.join(work, "uniform")
    make(p, False)
    d = max(np.abs(material_T(p, t) / T_COLD - 1).max() for t in (0.5, 1.0, 1.5))
    report("uniform temperature stays uniform", d < 1e-12, f"largest relative change {d:.1e}")
    log = ca.energy_log(p)
    want = -C_SAND * T_COLD * log["dis"][-1]
    r = abs(log["carried"][-1] / want - 1)
    report("heat carried out equals c T times the discharged mass", r < 1e-12,
           f"{log['carried'][-1]:.10e} J against {want:.10e} J ({r:.1e})")


def main(work):
    os.makedirs(work, exist_ok=True)
    t0 = time.perf_counter()
    check(work)
    print(f"STL heat checks took {time.perf_counter() - t0:.0f} s")
    if failed:
        print("FAILED: " + ", ".join(failed))
        sys.exit(1)
    print("test_stl_heat: PASS")


if __name__ == "__main__":
    main(sys.argv[1])
