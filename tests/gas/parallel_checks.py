"""WP F.9a: parallel runs with gas (decision 0134).

Why parallel runs with gas drift from serial, and the check that they are
right. Cases from tests/gas/fluid_cases.py (fluidised bed: sand d 0.2 mm,
8 x 40 cells, air ramped through the floor, model gasFlow with the pore gas
momentum, ergunWenYu) and tests/gas/hopper_gas_cases.py (plane-flow hopper,
d 0.1 mm, outlet opens at 0.15 s).

H  Processor boundaries: gasProperties parallelCheck compares, on every
   processor cut, the stored neighbour values of 15 cell fields (gas: p,
   Ug, rAU, H/A, e, k, phi_s, phi_g, void flag, q, inertia reaction; solid:
   alpha, rho, U, sigma) with freshly exchanged ones (a stale halo), and
   the two sides of 12 face fields (gas mass flux, face pressure, 1/A,
   H/A flux, old fluxes, inertia terms; solid: phiL, remap mass) with each
   other (a one-sided face value). Bed on 2 x 2 and 1 x 4 (cuts through the
   bed, the free surface and the void) to 0.03 s, every 5 steps; hopper
   (model gasFlow) on 2 x 2 to 0.2 s (outlet open), every 20 steps. Every
   mismatch within 1e-12 of the field's largest value.
A  Agreement before lift: the bed on 1, 2 x 2 and 1 x 4 processors to
   0.03 s. Gas pressure within 1e-10 of its largest |p|, grain velocity
   within 1e-8 of its largest |U| (the bed is nearly still: U is 1e-4 m/s
   and differs by 1e-13 m/s). Reported beside them, the round-off band: the
   same serial run on a renumbered mesh (renumberMesh, cells compared by
   their centres), which only changes the order of the sums. The gas
   density takes the absolute pressure p0 + p, so each step rounds the
   pressure by about eps p0 = 2e-11 Pa; solving the pressure or the void
   gas velocity to 1e-14 instead of 1e-12 changes nothing.
L  (--lift, tests/manual/test_lift_sensitivity.sh) After lift: serial
   against serial with only pTolerance changed (1e-12 to 1e-13 Pa), and
   serial against 2 x 1 and 1 x 4, to 0.28 s. Reported: the time where each
   pair first differs by 1e-3 of U's largest value. The serial pair parts
   too, so the parting after lift is the bed's own amplification of
   round-off, not a parallel error.

    python3 parallel_checks.py WORKDIR [--lift]
"""

import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import fluid_cases as fl  # noqa: E402
import gas_flow_checks as gf  # noqa: E402
import hopper_gas_cases as hc  # noqa: E402
import pore_cases as pc  # noqa: E402

NX = 8
N = NX * (fl.NB + fl.NV)
TOL = "pTolerance 1e-12;\n"
CHECK = "parallelCheck yes;\nparallelCheckEvery {};\n"
HALO_TOL = 1e-12
U_TOL = 1e-8
P_TOL = 1e-10


def worst_check(path):
    """Largest processor-boundary mismatch of the run, and the number of
    checks."""
    txt = open(f"{path}/log.celFoam").read()
    runs = re.findall(r"Parallel check \([^)]*\): largest \S+ \(\S+\), run (\S+)", txt)
    if not runs:
        return float("inf"), 0
    return float(runs[-1]), len(runs)


def bed(path, n, dec, *, end, write, extra=""):
    fl.make_bed(path, end=end, write_interval=write, decomp=dec, nx=NX, extra_gas=TOL + extra)
    return pc.run(path, n)


def centres_order(path):
    """Cell order sorted by centre (y, then x)."""
    C = gf.vfield(path, "0", "C", N)
    return np.lexsort((C[:, 0], C[:, 1]))


def compare(a, b, renumbered=False):
    """U and p over their largest values in the window. With renumbered,
    cells are matched by their centres."""
    oa = centres_order(a) if renumbered else slice(None)
    ob = centres_order(b) if renumbered else slice(None)
    du = dp = umax = pmax = 0.0
    for t in pc.times(a):
        ua, ub = gf.vfield(a, t, "U", N)[oa], gf.vfield(b, t, "U", N)[ob]
        pa, pb = pc.field(a, t, "p", N)[oa], pc.field(b, t, "p", N)[ob]
        umax = max(umax, float(np.abs(ua).max()))
        pmax = max(pmax, float(np.abs(pa).max()))
        du = max(du, float(np.abs(ua - ub).max()))
        dp = max(dp, float(np.abs(pa - pb).max()))
    return du / max(umax, 1e-300), dp / max(pmax, 1e-300)


def first_parting(a, b, level=1e-3):
    """First write time where U differs by level of its largest value
    over the run (None if never)."""
    ts = pc.times(a)
    umax = max(float(np.abs(gf.vfield(a, t, "U", N)).max()) for t in ts)
    for t in ts:
        d = float(np.abs(gf.vfield(a, t, "U", N) - gf.vfield(b, t, "U", N)).max())
        if d > level * umax:
            return float(t)
    return None


def main():
    work = sys.argv[1]
    lift = "--lift" in sys.argv[2:]
    os.makedirs(work, exist_ok=True)
    ok = True

    def report(name, good, text):
        nonlocal ok
        ok = ok and good
        print(f"{name}: {text}: {'PASS' if good else 'FAIL'}", flush=True)

    if lift:
        end, write = fl.T_RAMP, 0.005
        runs = {"serial": (1, (1, 1, 1), ""), "serial_tol13": (1, (1, 1, 1), "pTolerance 1e-13;\n"),
                "2x1": (2, (2, 1, 1), ""), "1x4": (4, (1, 4, 1), "")}
        for k, (n, dec, extra) in runs.items():
            bed(f"{work}/lift_{k}", n, dec, end=end, write=write, extra=extra)
        print("L  bed to 0.28 s (lift near 0.2 s): first write where U parts by 1e-3 of its largest value")
        for k in ("serial_tol13", "2x1", "1x4"):
            t = first_parting(f"{work}/lift_serial", f"{work}/lift_{k}")
            print(f"   serial vs {k}: {'never' if t is None else f'{t:.3f} s'}", flush=True)
        return 0

    # H: processor boundaries
    for k, (n, dec) in {"2x2": (4, (2, 2, 1)), "1x4": (4, (1, 4, 1))}.items():
        path = f"{work}/bed_{k}"
        bed(path, n, dec, end=0.03, write=0.005, extra=CHECK.format(5))
        w, nc = worst_check(path)
        report(f"H  bed {k}", w <= HALO_TOL and nc > 0,
               f"largest processor-boundary mismatch {w:.1e} over {nc} checks")
    path = f"{work}/hopper_2x2"
    hc.make_hopper(path, d=1e-4, end=0.2, write_interval=0.05, decomp=(2, 2, 1), rules=False,
                   extra_gas=TOL + CHECK.format(20))
    g = f"{path}/constant/gasProperties"
    pc.write(g, open(g).read().replace("model   poreGas;", "model   gasFlow;"))
    pc.run(path, 4, set_fields=True)
    w, nc = worst_check(path)
    report("H  hopper gasFlow 2x2", w <= HALO_TOL and nc > 0,
           f"largest processor-boundary mismatch {w:.1e} over {nc} checks")

    # A: agreement before lift
    for k in ("1", "renum"):
        path = f"{work}/bed_{k}"
        fl.make_bed(path, end=0.03, write_interval=0.005, nx=NX, extra_gas=TOL)
        pc.sh(path, ["blockMesh"], "blockMesh")
        if k == "renum":
            pc.sh(path, ["renumberMesh", "-overwrite"], "renumberMesh")
        pc.sh(path, ["postProcess", "-func", "writeCellCentres", "-time", "0"], "writeCellCentres")
        pc.sh(path, ["celFoam"], "celFoam")
    du, dp = compare(f"{work}/bed_1", f"{work}/bed_renum", renumbered=True)
    print(f"A  round-off band, serial vs serial renumbered: U {du:.1e}, p {dp:.1e} "
          "of their largest values", flush=True)
    for k in ("2x2", "1x4"):
        du, dp = compare(f"{work}/bed_1", f"{work}/bed_{k}")
        report(f"A  bed serial vs {k}", du <= U_TOL and dp <= P_TOL,
               f"U {du:.1e}, p {dp:.1e} of their largest values")
        drift = gf.drift(f"{work}/bed_{k}")
        report(f"   gas mass {k}", drift <= 1e-12, f"less the inflow constant to {drift:.1e}")

    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
