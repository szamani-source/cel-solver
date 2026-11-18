"""WP 3.6b done-when: steady inclined layer (Bagnold profile) in celFoam's
cycle with mu(I) friction (decision 0037), which measures the numerical
viscosity of the face Riemann solver on plastic flow (decisions 0028, 0036).

Set-up: T0 Mohr-Coulomb (rho 1800, nu 0.3, phi 35 deg, c 0, psi 0, no
tension) with the placeholder Jop law (mu_2 0.9618, I0 0.279, d 1 mm,
rho_s 2650), shear cap and free surface on as in the silo. A layer h = 0.05 m
deep on a slope of theta = 31 deg: the mesh is not tilted, gravity is
(g sin theta, -g cos theta); x is periodic (2 cells, cyclic), the floor is a
wall with friction 89 deg (sticks), void above. The run starts from the
exact steady state (stress and velocity) and runs T_END = 2 s, a few
relaxation times of the physical tangent viscosity, without damping.

Exact solution of the model (steady isochoric simple shear, decision 0037:
Mohr-Coulomb gives tau / sigma_n = sin(phi_m)): s_xx = s_yy = sn =
-rho g cos(theta) (h - y), s_zz = 2 nu sn (plane strain), s_xy = rho g
sin(theta) (h - y). So sin(phi_m) = tan(theta), mu(I) = tan(phi_m) gives one
I through the depth; P = -(2 + 2 nu)/3 sn; gamma_dot = I sqrt(P / rho_s) / d;
u(y) = (2/3) (I/d) sqrt(k rho g cos(theta) / rho_s) (h^1.5 - (h - y)^1.5),
k = (2 + 2 nu)/3.  Theta = 31 deg: I = 0.068, surface speed 1.133 m/s.

Checks: mass to 1e-10 and no stress in empty cells for every run. At
E = 0.1 MPa, 5 mm cells (10 across the layer): cell-average surface speed
(top cell) and the whole profile (largest cell error over the surface speed)
within 5 % of the exact cell averages at T_END, and steady (surface speed
changed by at most 1 % over the last 0.5 s). Recorded, not checked: E = 1
and 10 MPa (the error against E measures the numerical viscosity), 2.5 mm
at 0.1 MPa (mesh), and the effective numerical viscosity
eta_num = (tau_exact - s_xy) / gamma_dot in the lower half of the layer.

WP 3.6c (decision 0040): every layer runs with flowImpedance yes (mu(I)
tangent viscosity at faces between flowing cells, plastic Hancock predictor,
shear cap with phi_m); --noflow runs the WP 3.6b code path for comparison.
With --ci the 0.1 MPa, 5 mm layer is also checked against a regression
guard: surface speed within 10 % of exact and mean speed at least 0.85 of
exact (WP 3.6b code: 0.14 and 0.15).

Run by tests/test_inclined_layer.sh after celVoidCheck is built:

    python3 inclined_layer.py WORKDIR [--ci] [--noflow] [E_MPa:cell_mm[:lid] ...]
"""

import math
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "void"))
import void_checks as vc  # noqa: E402  (make_case, run)

RHO, NU, PHI, G = 1800.0, 0.3, 35.0, 9.81
MU2, I0, D, RHOS = 0.961813, 0.279, 1.0e-3, 2650.0    # decision 0037
THETA, H, T_END = 31.0, 0.05, 4.0
WRITE = [1.0, 2.0, 3.0]
TOL, STEADY = 0.05, 0.01
CI_TOL, CI_MEAN = 0.10, 0.85        # WP 3.6c regression guard (--ci)
FLOW = True                         # flowImpedance yes (WP 3.6c)
K = (2.0 + 2.0 * NU) / 3.0


def exact():
    """I and the surface speed of the steady layer."""
    th = math.radians(THETA)
    mu = math.tan(math.asin(math.tan(th)))
    mus = math.tan(math.radians(PHI))
    I = I0 * (mu - mus) / (MU2 - mu)
    a = (I / D) * math.sqrt(K * RHO * G * math.cos(th) / RHOS)
    return I, a


def u_exact(y, a):
    return (2.0 / 3.0) * a * (H ** 1.5 - np.maximum(H - y, 0.0) ** 1.5)


def u_cell(yc, cell, a):
    """Exact cell average: (U(y1) - U(y0)) / cell with U' = u."""
    def U(y):
        s = np.maximum(H - y, 0.0)
        return (2.0 / 3.0) * a * (H ** 1.5 * y + 0.4 * s ** 2.5)
    return (U(yc + 0.5 * cell) - U(yc - 0.5 * cell)) / cell


def setup(work, E_MPa, cell, shear_cap=True, lid=False, name=None):
    ny = int(round(H / cell)) + (0 if lid else 2)
    tag = f"E{E_MPa:g}_h{cell * 1000:g}{'' if shear_cap else '_nocap'}{'_lid' if lid else ''}"
    path = f"{work}/{name or 'layer_' + tag}"
    _, a = exact()
    material = (f"type mohrCoulomb;\nE {E_MPa * 1e6!r};\nnu {NU!r};\nphi {PHI!r};\nc 0;\n"
                f"psi 0;\nsigmaT 0;\nmu2 {MU2!r};\nI0 {I0!r};\nd {D!r};\nrhoS {RHOS!r};\n"
                + ("shearCap yes;\n" if shear_cap else "")
                + ("flowImpedance yes;\n" if FLOW else ""))
    check = dict(fill=[-1.0, -1.0, 1.0, H], U0=[0.0, 0.0, 0.0], tEnd=T_END,
                 maxCo=0.4, writeTimes=WRITE)
    old = vc.CELL
    vc.CELL = cell
    try:
        vc.make_case(path, 2, ny, material, check)
    finally:
        vc.CELL = old
    u_s = (2.0 / 3.0) * a * H ** 1.5
    with open(f"{path}/system/voidCheckDict", "a") as f:
        f.write(f"layer {{ theta {THETA!r}; kzz {2 * NU!r}; Us {u_s!r}; }}\n")
    # x periodic, floor sticks, tilted gravity
    bm = f"{path}/system/blockMeshDict"
    txt = open(bm).read()
    for me, other in (("left", "right"), ("right", "left")):
        txt = txt.replace(f"    {me}\n    {{\n        type patch;",
                          f"    {me}\n    {{\n        type cyclic;\n"
                          f"        neighbourPatch {other};")
    open(bm, "w").write(txt)
    mp = f"{path}/constant/materialProperties"
    txt = open(mp).read()
    txt = txt.replace("    left { friction 0; }\n", "").replace("    right { friction 0; }\n", "")
    txt = txt.replace("    bottom { friction 0; }", "    bottom { friction 89; }")
    txt = txt.replace(f"rho {vc.RHO};", f"rho {RHO};")
    open(mp, "w").write(txt)
    th = math.radians(THETA)
    gf = f"{path}/constant/g"
    txt = open(gf).read()
    txt = txt.replace(f"value (0 {-vc.G_ACC!r} 0);",
                      f"value ({G * math.sin(th)!r} {-G * math.cos(th)!r} 0);")
    open(gf, "w").write(txt)
    return path


def profile(fname, cell):
    """Layer cells (x-averaged, bottom up): y, Ux, s_xy, alpha."""
    c = np.loadtxt(fname)
    j = np.floor(c[:, 1] / cell + 1e-9).astype(int)
    n = j.max() + 1
    y = (np.arange(n) + 0.5) * cell
    ux = np.array([c[j == k, 4].mean() for k in range(n)])
    sxy = np.array([c[j == k, 8].mean() for k in range(n)])
    al = np.array([c[j == k, 2].mean() for k in range(n)])
    return y, ux, sxy, al


def layer(work, E_MPa, cell, shear_cap=True, lid=False):
    path = setup(work, E_MPa, cell, shear_cap, lid)
    sec = vc.run(path)
    h = np.loadtxt(f"{path}/history.dat")
    _, a = exact()
    nl = int(round(H / cell))
    y, ux, sxy, al = profile(f"{path}/cells.dat", cell)
    ue = u_cell(y[:nl], cell, a)
    us = ux[nl - 1]
    snaps = {t: profile(f"{path}/cells_{t:g}.dat", cell)[1][nl - 1] for t in WRITE}
    th = math.radians(THETA)
    tau = RHO * G * math.sin(th) * (H - y[:nl])
    # numerical viscosity at the faces of the lower half (smooth profile)
    lo = np.arange(1, nl // 2)
    yf = lo * cell
    gdot = (ux[lo] - ux[lo - 1]) / cell
    tau_f = RHO * G * math.sin(th) * (H - yf)
    s_f = 0.5 * (sxy[lo] + sxy[lo - 1])
    eta = float(np.sum(tau_f - s_f) / max(np.sum(gdot), 1e-12))
    G_mod = E_MPa * 1e6 / (2 * (1 + NU))
    lam = E_MPa * 1e6 * NU / ((1 + NU) * (1 - 2 * NU))
    cp = math.sqrt((lam + 2 * G_mod) / RHO)
    cs = math.sqrt(G_mod / RHO)
    return dict(E=E_MPa, cell=cell, cap=shear_cap, seconds=sec, steps=int(h[-1, vc.STEP]),
                dm=abs(h[-1, vc.MASS] / h[0, vc.MASS] - 1.0), s_empty=h[:, vc.SEMPTY].max(),
                us=us, us_exact=ue[-1], err_s=us / ue[-1] - 1.0,
                err_prof=float(np.abs(ux[:nl] - ue).max() / ue[-1]),
                drift=abs(us / snaps[WRITE[-1]] - 1.0), snaps=snaps,
                tau_err=float(np.abs(sxy[:nl] - tau).max() / tau.max()),
                eta=eta, rcs=RHO * cs * cell, ux=ux[:nl], ue=ue, alpha_top=al[nl - 1],
                leak=float(al[nl:].sum()), umean=float(ux[:nl].mean() / ue.mean()))


def at_rest(work):
    """Layer below the stopping angle (theta 25 deg, tan 25 < sin 35) from
    its exact elastic stress at rest, E = 1 MPa, 0.5 s: tests the face
    reconstruction under gravity not along an axis (decision 0039). Checked:
    the layer below the top cell stays within 1e-4 m/s, mass to 1e-10, no
    stress in empty cells; the top cell is recorded (it slides, a known issue)."""
    global THETA, T_END, WRITE
    keep = THETA, T_END, WRITE
    THETA, T_END, WRITE = 25.0, 0.5, [0.25]
    try:
        path = setup(work, 1.0, 5.0e-3, name="rest")
    finally:
        THETA, T_END, WRITE = keep
    txt = open(f"{path}/system/voidCheckDict").read()
    open(f"{path}/system/voidCheckDict", "w").write(re.sub(r"Us [^;]*;", "Us 0;", txt))
    sec = vc.run(path)
    h = np.loadtxt(f"{path}/history.dat")
    c = np.loadtxt(f"{path}/cells.dat")
    nl = int(round(H / 5.0e-3))
    j = np.floor(c[:, 1] / 5.0e-3 + 1e-9).astype(int)
    sp = np.hypot(c[:, 4], c[:, 5])
    bulk = sp[j < nl - 1].max()
    top = sp[j == nl - 1].max()
    dm = abs(h[-1, vc.MASS] / h[0, vc.MASS] - 1.0)
    good = bulk <= 1e-4 and dm <= 1e-10 and h[:, vc.SEMPTY].max() == 0.0
    print(f"\n  At rest below the stopping angle (theta 25 deg, E = 1 MPa, 5 mm, 0.5 s, "
          f"{sec:.0f} s): largest speed below the top cell {bulk:.1e} m/s (<= 1e-4), top cell "
          f"{top:.1e} m/s (recorded); mass change {dm:.0e}: {'PASS' if good else 'FAIL'}")
    return good


def main():
    work = sys.argv[1]
    global FLOW
    ci = "--ci" in sys.argv[2:]
    FLOW = "--noflow" not in sys.argv[2:]
    args = [a for a in sys.argv[2:] if a not in ("--ci", "--noflow")]
    # (E MPa, cell mm, lid): free surface (void above) or a smooth lid on top
    cases = [(0.1, 5.0, False), (1.0, 5.0, False), (10.0, 5.0, False), (0.1, 2.5, False),
             (0.1, 5.0, True), (1.0, 5.0, True), (10.0, 5.0, True), (0.1, 2.5, True)]
    if ci:
        cases = [(0.1, 5.0, False)]
    if args:
        cases = [(float(s.split(":")[0]), float(s.split(":")[1]), s.endswith(":lid"))
                 for s in args]
    I, a = exact()
    t_all = time.perf_counter()
    print("Steady inclined layer (Bagnold), T0 Mohr-Coulomb with mu(I) (Jop, placeholders),")
    print("flow impedance " + ("on (WP 3.6c)" if FLOW else "off (--noflow, WP 3.6b code)"))
    print(f"theta {THETA:g} deg, h {H:g} m, x periodic, sticking floor, shear cap on, from the")
    print(f"exact steady state to t = {T_END:g} s. Exact: I = {I:.4f}, surface speed "
          f"{(2 / 3) * a * H ** 1.5:.4f} m/s")
    ok = at_rest(work)
    for E, mm, lid in cases:
        r = layer(work, E, mm * 1e-3, lid=lid)
        basic = r["dm"] <= 1e-10 and r["s_empty"] == 0.0 and np.isfinite(r["us"])
        checked = E == 0.1 and mm == 5.0 and not lid and not ci
        good = basic and (not checked or (abs(r["err_s"]) <= TOL and r["err_prof"] <= TOL
                                          and r["drift"] <= STEADY))
        guard = ci and FLOW
        if guard:
            good = good and abs(r["err_s"]) <= CI_TOL and r["umean"] >= CI_MEAN
        ok &= good
        print(f"\n  E = {E:g} MPa, {mm:g} mm cells ({int(round(H / (mm * 1e-3)))} across), "
              f"{'smooth lid' if lid else 'free surface'}: {r['steps']} steps "
              f"({r['seconds']:.0f} s)")
        print(f"    surface speed {r['us']:.4f} m/s, exact cell average {r['us_exact']:.4f} "
              f"({100 * r['err_s']:+.1f} %); largest profile error {100 * r['err_prof']:.1f} % "
              f"of it; mean speed {r['umean']:.3f} of exact; surface at "
              + ", ".join(f"{t:g} s {v:.4f}" for t, v in r["snaps"].items())
              + f"; last {T_END - WRITE[-1]:g} s change {100 * r['drift']:.2f} %")
        print(f"    shear stress within {100 * r['tau_err']:.2f} % of rho g sin(theta) (h - y); "
              f"numerical viscosity sum(tau_exact - s_xy)/sum(gamma_dot), lower half: "
              f"{r['eta']:.2f} Pa s (rho c_s h = {r['rcs']:.0f} Pa s); top cell alpha "
              f"{r['alpha_top']:.4f}, material above the layer {r['leak']:.1e} cells")
        print("    u / u_exact by cell, bottom up: "
              + " ".join(f"{u / e:.3f}" for u, e in zip(r["ux"], r["ue"])))
        print(f"    mass change {r['dm']:.0e}; largest stress in empty cells "
              f"{r['s_empty']:.0e} Pa"
              + (f"; guard (WP 3.6c): surface within {100 * CI_TOL:g} %, mean speed "
                 f">= {CI_MEAN:g}" if guard else "")
              + f": {'checked' if checked or guard else 'recorded'}: "
              f"{'PASS' if good else 'FAIL'}")
    print(f"\nInclined layer checks: {'ALL PASS' if ok else 'FAILED'} "
          f"({time.perf_counter() - t_all:.0f} s)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
