"""WP 3.8b done-when: the free surface (freeSurface yes, decision 0031) on an
axisymmetric wedge mesh (decision 0061), as tests/test_void.sh on the
Cartesian mesh.

Wedge 2 x 2.5 deg about the y axis, radial direction x, 10 mm cells, axis
cells are prisms. Material T0 Mohr-Coulomb (rho 1800, E 10 MPa, nu 0.3,
phi 35) unless said.

1. Lake at rest: K0 column between a smooth side wall and floor, R = 0.1 m,
   in a 0.6 m tall cylinder filled to 0.455 m (top cell half full), full
   gravity, no damping, 1000 steps: speeds <= 1e-8 m/s, no stress in empty
   cells, mass to 1e-12.
2. Falling plug: an elastic cylinder r <= 0.05 m, 0.1 m tall, falls 0.1 m
   through void in a cylinder R = 0.1 m onto a rigid floor: acceleration
   within 0.5 % of g and |stress| <= 1e-6 Pa while falling; lowest centroid
   >= 0.04 m, peak floor force above the weight, mass to 1e-12.
3. Void-top tank: cylindrical tank R = 0.2 m, 20 deg wall and floor, filled
   to H = 1.2 m with void above, T0 rule set, damped from zero stress to
   0.5 s (ramp, damping as the WP 3.8a Janssen tank), then 10 % more steps
   undamped: stress change <= 2 % of rho g H and KE / strain energy <= 1e-5;
   mean s_v below depth 2R within 5 % of circular Janssen with K measured
   at the wall (K in 0.27 to 0.45); mass to 1e-12.
4. Two processors (split across the radius) match serial to 1e-10 of
   rho g H over 300 steps of case 3 without the T0 rules (the normal cap's
   bisection is not bit-for-bit across processors, a known issue of WP 3.4b). (The elastic plug's landing is not
   used: waves focusing on the axis make it chaotic, round-off grows to
   1e-3 after the impact.)

Run by tests/test_wedge_void.sh after celVoidCheck and celSiloCheck are
built:

    python3 wedge_void_checks.py WORKDIR
"""

import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
for d in ("void", "silo", "lagrange", "wedge"):
    sys.path.insert(0, os.path.join(HERE, "..", d))
import cases as cs  # noqa: E402
import void_checks as vc  # noqa: E402
import silo_checks as sc  # noqa: E402
import wedge_checks as wc  # noqa: E402

CELL = 0.01
RHO, G_ACC = vc.RHO, vc.G_ACC
STEP, T, MASS, PX, PY, CX, CY, SMAX, SEMPTY, UMAX, FY = range(11)
COS_T = math.cos(math.radians(2.5))   # wedge vertices sit at x = r cos(2.5 deg)


def wedge_void_case(path, R, H, material, check, gravity=True, nProcs=1):
    p = vc.make_case(path, 1, 1, material, check, gravity, nProcs)
    wc.wedge_mesh(p, 0.0, R, H, CELL)
    mp = f"{p}/constant/materialProperties"
    txt = open(mp).read()
    txt = txt.replace("    left { friction 0; }\n", "")
    open(mp, "w").write(txt)
    return p


def check_lake(work):
    p = wedge_void_case(f"{work}/lake", 0.1, 0.6, vc.MC,
                        {"fill": [-1.0, -1.0, 1.0, 0.455], "U0": [0, 0, 0],
                         "initK0": vc.NU / (1 - vc.NU), "nSteps": 1000})
    t = vc.run(p)
    h = vc.history(p)
    vmax = h[1:, UMAX].max()
    dm = abs(h[-1, MASS] / h[0, MASS] - 1)
    se = h[:, SEMPTY].max()
    ok = vmax <= 1e-8 and se == 0.0 and dm <= 1e-12
    return sc.report("Lake at rest in a cylinder (R 0.1 m, filled to 0.455 m, top cell "
                     "half full), 1000 steps", ok,
                     [f"largest speed {vmax:.1e} m/s (<= 1e-8); stress in empty cells "
                      f"{se:.1e} Pa (0); mass change {dm:.0e} (<= 1e-12); {t:.0f} s"])


def falling(work, name, nProcs=1):
    y0 = 0.10
    t_land = math.sqrt(2 * y0 / G_ACC)
    p = wedge_void_case(f"{work}/{name}", 0.1, 0.3, vc.ELASTIC,
                        {"fill": [-1.0, y0, 0.05 * COS_T, y0 + 0.1], "U0": [0, 0, 0],
                         "tEnd": t_land + 0.03}, nProcs=nProcs)
    t = vc.run(p, nProcs)
    return p, t, y0, t_land


def check_fall(work):
    p, t, y0, t_land = falling(work, "fall")
    h = vc.history(p)
    fall = (h[:, T] > 0) & (h[:, T] < 0.9 * t_land)
    acc = -2 * np.polyfit(h[fall, T], h[fall, CY], 2)[0]
    e_acc = abs(acc / G_ACC - 1)
    s_fall = h[fall, SMAX].max()
    after = h[:, T] >= 0.9 * t_land
    ymin = h[after, CY].min()
    weight = h[0, MASS] * G_ACC
    fb = h[after, FY].max()
    dm = abs(h[-1, MASS] / h[0, MASS] - 1)
    ok = e_acc <= 0.005 and s_fall <= 1e-6 and ymin >= 0.04 and fb > weight and dm <= 1e-12
    return sc.report("Elastic plug (r 0.05 m, 0.1 m tall) falling 0.1 m onto the floor", ok,
                     [f"centroid acceleration {acc:.5f} m/s^2, {100 * e_acc:.3f} % from g "
                      f"(<= 0.5 %); largest |stress| while falling {s_fall:.1e} Pa (<= 1e-6)",
                      f"landing at {math.sqrt(2 * G_ACC * y0):.2f} m/s: lowest centroid "
                      f"{ymin:.4f} m (>= 0.04); peak floor force {fb / weight:.1f} times the "
                      f"weight (> 1); mass change {dm:.0e} (<= 1e-12); {t:.0f} s"])


TANK_R, TANK_H, TANK_HM = 0.2, 1.2, 1.25
TANK_T = 0.5        # damped to this time, as the WP 4.3c void-top tank
# T0 rule set (decision 0050): mu(I) placeholders, shear cap, rule 9, flow
# impedance (void below 1e-4 counts as full), normal cap. Without it, void
# round-off (1e-9) in a full region makes partly filled cells that never
# settle, on the Cartesian mesh too
T0_RULES = ("mu2 0.961813;\nI0 0.279;\nd 0.001;\nrhoS 2650;\nshearCap yes;\n"
            "spreadLayers yes;\nflowImpedance yes;\nnormalCap yes;\n")


def tank_case(work, name, check, nProcs=1, rules=True):
    """Cylindrical tank (wedge), filled to TANK_H with void above; with
    nProcs 2 the split is across the radius."""
    p = sc.make_case(f"{work}/{name}", TANK_R, TANK_HM, True, check, sc.T_RAMP, sc.ALPHA,
                     nProcs)
    wc.wedge_mesh(p, 0.0, TANK_R, TANK_HM, CELL)
    wc.set_walls(p, {"right": wc.DELTA_W, "bottom": wc.DELTA_W, "top": 0.0})
    open(f"{p}/constant/materialProperties", "a").write("freeSurface yes;\n" + (T0_RULES if rules else ""))
    open(f"{p}/system/siloCheckDict", "a").write(f"fillBelow (0 1 0 {TANK_H!r});\n")
    cs.write(f"{p}/system/decomposeParDict",
             cs.HEADER.format(cls="dictionary", obj="decomposeParDict")
             + f"numberOfSubdomains {nProcs};\nmethod simple;\ncoeffs {{ n ({nProcs} 1 1); }}\n")
    return p


def check_tank(work):
    R, H, Hm = TANK_R, TANK_H, TANK_HM
    p = tank_case(work, "tank", dict(mode="settle", keTol=sc.KE_TOL, checkEvery=20, tMax=TANK_T,
                                     continueFraction=0.1))
    sec = sc.run(p)
    s = sc.summary(p)
    _, det = sc.settle_report(s, RHO * G_ACC * H)
    # The surface cells keep creeping (about 1 mm/s), so the KE criterion of
    # settle mode is not met by TANK_T; the undamped continuation must hold
    # the stresses, as the Cartesian void-top tank (WP 4.3c)
    good = (s["contMaxStressChange"] / (RHO * G_ACC * H) <= 0.02
            and s["contKEratio"] <= 10 * sc.KE_TOL)
    mu = math.tan(math.radians(wc.DELTA_W))
    c = sc.cells(p)
    nr, ny = int(round(R / CELL)), int(round(Hm / CELL))
    g = {k: c[:, i].reshape(ny, nr) for i, k in
         enumerate(("x", "y", "ux", "uy", "sxx", "syy", "szz", "sxy", "rho", "ep"))}
    z = H - g["y"][:, 0]
    r = g["x"][0]
    w = r / r.sum()
    sv = -(g["syy"] * w).sum(axis=1)
    deep = z >= 2 * R
    sh = -(1.5 * g["sxx"][:, -1] - 0.5 * g["sxx"][:, -2])
    K = sh[deep].sum() / sv[deep].sum()
    err = np.abs(wc.janssen_circ(z, R, K, mu)[deep] / sv[deep] - 1.0).max()
    dm = abs(s["massChange"])
    ok = good and err <= 0.05 and 0.27 <= K <= 0.45 and dm <= 1e-12
    det += [f"K at the wall {K:.4f} (0.27 to 0.45); mean s_v against circular Janssen, "
            f"depth >= 2R: largest error {100 * err:.2f} % (<= 5 %)"]
    for frac in (0.25, 0.5, 1.0):
        i = int(np.argmin(np.abs(z - frac * H)))
        det.append(f"  depth {z[i]:.3f} m: s_v {sv[i]:8.1f} Pa, Janssen "
                   f"{wc.janssen_circ(z[i], R, K, mu):8.1f}")
    det.append(f"mass change {dm:.1e} (<= 1e-12); wall clock {sec:.0f} s")
    return sc.report("Cylindrical tank with a void top (R 0.2 m, filled to 1.2 m, 20 deg)",
                     ok, det)


def check_parallel(work):
    chk = dict(mode="fixed", nSteps=300)
    a = tank_case(work, "par1", chk, rules=False)
    b = tank_case(work, "par2", chk, nProcs=2, rules=False)
    sc.run(a)
    t2 = sc.run(b, 2)
    ca, cb = sc.cells(a), sc.cells(b, 2)
    d = np.abs(ca[:, 4:8] - cb[:, 4:8]).max() / (RHO * G_ACC * TANK_H)
    return sc.report("Two processors against serial (void-top tank, 300 steps, split across "
                     "the radius)", d <= 1e-10,
                     [f"largest stress difference {d:.1e} of rho g H (<= 1e-10); {t2:.0f} s"])


def main():
    work = sys.argv[1]
    t0 = time.perf_counter()
    res = [check_lake(work), check_fall(work), check_tank(work), check_parallel(work)]
    ok = all(res)
    print(f"\nWedge free-surface checks: {'ALL PASS' if ok else 'FAILED'} "
          f"({time.perf_counter() - t0:.0f} s)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
