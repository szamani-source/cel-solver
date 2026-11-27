"""WP F.8i: the pore pressure drains at the material's own surface
(gasProperties frictionOnly with diffusivity; decision 0133).

A still, fully supported layer of the dam-break grains (tests/gas/
collapse_cases.py) lies on a sealed floor under gravity, 0.1 m wide, its
top cell partly filled. With frictionOnly the pressure diffuses with the
given kappa and is zero at the surface, so the pressure in the bottom cell
follows the series for a layer of depth h sealed at the base:

    p(z, t) / (rho_b g h) = sum 2 (-1)^n / (h beta_n)^2 sin(beta_n z) exp(-kappa beta_n^2 t),
    beta_n = (2n + 1) pi / (2 h), z the depth.

D  layers whose top cell is 1/2, 3/4 and 1/4 full (10 mm cells, h 25, 27.5
   and 22.5 mm) and 1/2 full on 5 mm cells (h 22.5 mm): the bottom-cell
   pressure within 2 % of rho_b g h of the series at t = h^2/kappa times 0.25,
   0.5 and 1. (Before WP F.8i the zero pressure sat on the face above the
   top cell: 15 % to 28 % too slow with a half-filled top cell.)
P  the 10 mm, h 25 mm layer on 2 processors (cut across the columns): p
   within 1e-10 of rho_b g h of serial to t = h^2/(2 kappa). Later the solid
   itself parts serial from parallel by round-off (a cohesionless layer at the
   Mohr-Coulomb apex, the known F.9 issue), so those writes are reported.

    python3 surface_drainage_checks.py WORKDIR
"""

import math
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import collapse_cases as cc  # noqa: E402

W = 0.1                 # layer width [m]
KAPPA = 0.11            # about k p0 / (mu e) of the dam-break grains [m^2/s]
res = []


def report(name, good, detail):
    res.append(bool(good))
    print(f"  {name}: {detail}: {'PASS' if good else 'FAIL'}", flush=True)


def setup(cell, h):
    cc.CELL = cell
    cc.NX = int(round(W / cell))
    cc.NY = int(math.ceil((h + 0.03) / cell - 1e-9))
    cc.N = cc.NX * cc.NY
    cc.X0, cc.H0 = W, h
    cc.supply = lambda: 0.0                 # sealed floor
    cc.WRITE = round(h * h / KAPPA / 4, 9)


def make(path, cell, h, nProcs=1):
    setup(cell, h)
    cc.make(path, decomp=(nProcs, 1, 1), friction_only=True, kappa=KAPPA, end=4 * cc.WRITE,
            extra_gas="pTolerance 1e-12;\n")
    # partly filled top cell, isotropic lithostatic total stress, p = weight
    alpha, rho, sig, p = [], [], [], []
    for j in range(cc.NY):
        yb = j * cell
        for i in range(cc.NX):
            f = min(max((h - yb) / cell, 0.0), 1.0)
            w = cc.RHO_B * cc.G_ACC * ((h - yb - 0.5 * cell) if f >= 1 else 0.5 * f * cell)
            alpha.append(f)
            rho.append(cc.RHO_B * f)
            sig.append((-w, -w))
            p.append(w if f > 0 else 0.0)
    cc.fl._write_fields(path, alpha, rho, sig, cc.N)
    f = f"{path}/0/p"
    txt = re.sub(r"internalField nonuniform List<scalar> \d+\n\(.*?\n\)",
                 lambda m: f"internalField nonuniform List<scalar> {cc.N}\n(\n"
                           + "\n".join(repr(v) for v in p) + "\n)", open(f).read(), flags=re.S)
    cc.pc.write(f, txt)
    # half of the floor stays the sealed floorRes patch, so both patches exist
    f = f"{path}/system/topoSetDict"
    cc.pc.write(f, open(f).read().replace(f"({W!r} 0.001 1)", f"({W / 2!r} 0.001 1)"))


def series(h, z, t):
    s = 0.0
    for n in range(400):
        b = (2 * n + 1) * math.pi / (2 * h)
        s += 2 * (-1) ** n / (h * b) ** 2 * math.sin(b * z) * math.exp(-KAPPA * b * b * t)
    return s


def writes(path):
    return sorted((d for d in os.listdir(path) if re.fullmatch(r"[0-9.e-]+", d) and float(d) > 0), key=float)


def main():
    work = sys.argv[1]
    t_all = time.perf_counter()
    print(f"D  still layer {W} m wide, sealed floor, frictionOnly, kappa {KAPPA} m^2/s; "
          f"bottom-cell p against the series")
    for cell, h, top in ((0.01, 0.025, "1/2"), (0.01, 0.0275, "3/4"), (0.01, 0.0225, "1/4"),
                         (0.005, 0.0225, "1/2")):
        path = f"{work}/D_{cell}_{h}"
        make(path, cell, h)
        secs = cc.run(path, 1)
        w0 = cc.RHO_B * cc.G_ACC * h
        worst, parts = 0.0, []
        for t in writes(path):
            p = cc.pc.field(path, t, "p", cc.N).reshape(cc.NY, cc.NX)[0, cc.NX // 2] / w0
            ref = series(h, h - 0.5 * cell, float(t))
            worst = max(worst, abs(p - ref))
            parts.append(f"{p:.3f}/{ref:.3f}")
        report(f"D {cell * 1000:g} mm cells, h {h * 1000:g} mm (top cell {top} full)", worst <= 0.02,
               f"p/(rho_b g h) at t/(h^2/kappa) 0.25 to 1: {', '.join(parts)}; largest error "
               f"{100 * worst:.2f} % (<= 2 %); {secs:.0f} s")
    cell, h = 0.01, 0.025
    path = f"{work}/P_2"
    make(path, cell, h, nProcs=2)
    cc.run(path, 2)
    w0 = cc.RHO_B * cc.G_ACC * h
    diffs = [np.abs(cc.pc.field(path, t, "p", cc.N) - cc.pc.field(f"{work}/D_{cell}_{h}", t, "p", cc.N)).max() / w0
             for t in writes(path)]
    report("P 2 processors", max(diffs[:2]) <= 1e-10,
           f"p within {max(diffs[:2]):.1e} of rho_b g h of serial to h^2/(2 kappa) (<= 1e-10); later "
           + ", ".join(f"{d:.1e}" for d in diffs[2:]) + " (the unsupported layer is at the Mohr-Coulomb "
           "apex, where the solid parts by round-off, F.9)")
    print(f"\ntest_surface_drainage: {'PASS' if all(res) else 'FAIL'} ({sum(res)}/{len(res)}), "
          f"{time.perf_counter() - t_all:.0f} s")
    return 0 if all(res) else 1


if __name__ == "__main__":
    sys.exit(main())
