"""WP C.2 done-when checks: a cohesive column with a free vertical face
(decision 0112). About 4 minutes on 4 cores.

The a = 0.5 column of tests/collapse (H = 0.1 m, L = 0.2 m against the
left wall, rough floor, rho 1800, E 1 MPa, phi 35, psi 0, shear cap and
rule 9) with cohesion c and the default tensile strength c cot(phi), run
to 0.5 s. It stands if its kinetic energy never exceeds 1e-5 J/m (the
elastic settling of a standing column peaks at 1e-6 to 4e-6 J/m; a column
that fails passes 3e-5 J/m within 0.5 s); run-out is reported.

Reference: Terzaghi's plane-wedge critical height of a vertical cut,
H_c = n (c / gamma) tan(45 + phi/2), n = 4 without a tension crack and 2.67
with one, so the critical cohesion for H = 0.1 m lies between
c = gamma H / (4 tan(62.5)) = 230 Pa and gamma H / (2.67 tan(62.5)) = 344 Pa.

[1] Mohr-Coulomb on 5 and 2.5 mm: c = 250 Pa fails and 300 Pa stands on
    both, so the critical c is in (250, 300) Pa on both meshes: inside
    Terzaghi's limits and the same within 20 %; c = 275 Pa reported, which
    narrows the bracket to 10 % on each mesh. 10 mm cells reported.
[2] Drucker-Prager (plane-strain fit, as the silo, decision 0102) on 5 mm:
    c = 230 Pa fails and 344 Pa stands (critical c inside the limits).

    python3 vertical_cut.py WORKDIR
"""

import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "collapse"))
sys.path.insert(0, os.path.join(HERE, "..", "void"))
import collapse_checks as cc  # noqa: E402
import void_checks as vc  # noqa: E402

KE_STAND = 1e-5
T_END = 0.5
fails = []


def report(name, ok, detail):
    print(f"  {name}: {'PASS' if ok else 'FAIL'} ({detail})", flush=True)
    if not ok:
        fails.append(name)


def one(args):
    work, model, h, c = args
    vc.CELL = h
    cc.CELL = h
    L_i, H_i, nx, ny = cc.geometry(0.5)
    p = cc.setup(work, 0.5, 1.0, T_END, f"{model}_c{c:g}_{h * 1000:g}mm")
    mp = f"{p}/constant/materialProperties"
    txt = open(mp).read().replace("c 0;", f"c {c!r};").replace("sigmaT 0;\n", "")
    if model == "DP":
        txt = txt.replace("type mohrCoulomb;", "type druckerPrager;\nfit plane_strain;")
    open(mp, "w").write(txt)
    vc.run(p)
    hist = np.loadtxt(f"{p}/history.dat")
    r = cc.measure(f"{p}/cells.dat", nx, ny, L_i)
    return (model, h, c), (hist[:, cc.KE].max(), r["runout"])


def main():
    work = os.path.abspath(sys.argv[1])
    os.makedirs(work, exist_ok=True)
    t0 = time.perf_counter()
    jobs = [("MC", h, c) for h in (0.0025, 0.005) for c in (250.0, 275.0, 300.0)]
    jobs += [("MC", 0.01, c) for c in (250.0, 300.0)]
    jobs += [("DP", 0.005, c) for c in (230.0, 344.0)]
    with ProcessPoolExecutor(4) as ex:
        res = dict(ex.map(one, [(work, *j) for j in jobs]))
    stands = {k: v[0] <= KE_STAND for k, v in res.items()}
    for (m, h, c), (ke, ro) in sorted(res.items()):
        print(f"  {m} {h * 1000:4g} mm c {c:4.0f} Pa: KE peak {ke:.1e} J/m, run-out {ro:.3f}: "
              f"{'stands' if stands[(m, h, c)] else 'fails'}", flush=True)
    print("[1] Mohr-Coulomb, critical c in (250, 300) Pa on 5 and 2.5 mm "
          "(Terzaghi 230 to 344 Pa)", flush=True)
    for h in (0.005, 0.0025):
        lo, hi = not stands[("MC", h, 250.0)], stands[("MC", h, 300.0)]
        mid = "stands" if stands[("MC", h, 275.0)] else "fails"
        report(f"{h * 1000:g} mm", lo and hi,
               f"250 Pa {'fails' if lo else 'stands'}, 300 Pa {'stands' if hi else 'fails'}; "
               f"275 Pa {mid}")
    print("[2] Drucker-Prager (plane strain fit), 5 mm", flush=True)
    lo, hi = not stands[("DP", 0.005, 230.0)], stands[("DP", 0.005, 344.0)]
    report("Terzaghi limits", lo and hi,
           f"230 Pa {'fails' if lo else 'stands'}, 344 Pa {'stands' if hi else 'fails'}")
    print(f"{len(fails)} failures ({time.perf_counter() - t0:.0f} s)")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
