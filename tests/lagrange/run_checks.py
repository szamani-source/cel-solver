#!/usr/bin/env python3
"""WP 2.3 done-when checks for eulerianSolid::lagrangeStep (C++).

Run by tests/test_lagrange_step.sh after celLagrangeStepCheck is built:

    python3 run_checks.py WORKDIR [MAXCO]

MAXCO is the default maxCo of eulerianSolid (read by the shell script from
eulerianSolid.C), at which the stability cases must be stable.
Material: rho 2000 kg/m^3, E 20 MPa, nu 0.3 (decision 0013). Waves use
strain amplitude 1e-6 (no remap yet, decision 0014). Prints one PASS/FAIL
line per check and exits non-zero on any failure.
"""

import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "prototype", "constitutive"))

import cases as cs  # noqa: E402
import models  # noqa: E402

K, G, M = cs.moduli()
CP = math.sqrt(M / cs.RHO)
CS = math.sqrt(G / cs.RHO)
A = 1e-6

failures = []


def report(name, ok, detail):
    print(f"{'PASS' if ok else 'FAIL'}  {name}: {detail}", flush=True)
    if not ok:
        failures.append(name)


def polarisation(mode, k):
    kv = np.array(k, float)
    kh = kv / np.linalg.norm(kv)
    if mode == "P":
        return kh
    p = np.array([0.0, 1.0, 0.0]) if k[1] == 0 else np.array([-1.0, 1.0, 0.0])
    p -= p.dot(kh) * kh
    return p / np.linalg.norm(p)


def wave_case(work, name, n, mode="P", k=(1, 0, 0), profile="sine",
              crossings=1.0, dim=2, maxCo=None, secondOrder=True, nProcs=1,
              split=None, sides=None, ny=None, pulse=None):
    check = dict(setup="wave", mode=mode, k=list(k), profile=profile,
                 amplitude=A, L=1.0, n=n, crossings=crossings,
                 polarisation=list(polarisation(mode, k)))
    if pulse is not None:
        check["pulse"] = list(pulse)
    path = cs.make_case(os.path.join(work, name), n, sides or {}, check,
                        dim=dim, maxCo=maxCo, secondOrder=secondOrder,
                        nProcs=nProcs, split=split, ny=ny)
    return cs.run(path, nProcs)


def momentum_drift(path):
    h = cs.history(path)
    return np.abs(h[:, 10:13] - h[0, 10:13]).max() / h[:, 13].max()


def sine_error(path, mode, k):
    """L1 error of the polarised velocity against the exact sine wave, and
    the phase speed of the fundamental relative to the exact speed."""
    h = cs.history(path)
    c = cs.cells(path)
    t = h[-1, 1]
    kv = np.array(k, float)
    kn = np.linalg.norm(kv)
    cw = CP if mode == "P" else CS
    v = c[:, 3:6] @ polarisation(mode, k)
    phi = 2 * np.pi * (c[:, :3] @ kv)
    wt = 2 * np.pi * kn * cw * t
    exact = -cw * A * np.sin(phi - wt)
    err = np.abs(v - exact).mean() / (cw * A)
    # v = -c A sin(phi - w t): Fourier coefficient (i/2) c A exp(-i w t)
    a = np.sum(v * np.exp(-1j * phi))
    ph = -np.angle(a / (0.5j * cw * A))
    ph = wt + (ph - wt + np.pi) % (2 * np.pi) - np.pi
    speed = ph / (2 * np.pi * kn * t)
    return err, speed / cw - 1


def check_waves(work):
    """P and S waves along x and at 45 degrees, periodic: speed within 1 %,
    order >= 1.7 over 3 meshes, momentum and rho V_L = m."""
    worst_drift, worst_rho = 0.0, 0.0
    for mode, k in [("P", (1, 0, 0)), ("S", (1, 0, 0)),
                    ("P", (1, 1, 0)), ("S", (1, 1, 0))]:
        errs, speeds = [], []
        for n in (25, 50, 100):
            p = wave_case(work, f"wave{mode}{k[1]}_{n}", n, mode, k)
            e, s = sine_error(p, mode, k)
            errs.append(e)
            speeds.append(s)
            worst_drift = max(worst_drift, momentum_drift(p))
            worst_rho = max(worst_rho, cs.history(p)[:, -1].max())
        orders = [math.log2(errs[0] / errs[1]), math.log2(errs[1] / errs[2])]
        angle = "45 deg" if k[1] else "along x"
        report(f"{mode}-wave {angle} speed",
               max(abs(s) for s in speeds) <= 0.01,
               "speed error " + ", ".join(f"{100 * s:+.3f} %" for s in speeds)
               + " on 25, 50, 100 cells (<= 1 %)")
        report(f"{mode}-wave {angle} order",
               min(orders) >= 1.7,
               f"L1 errors {errs[0]:.3e} {errs[1]:.3e} {errs[2]:.3e}, orders "
               f"{orders[0]:.2f} {orders[1]:.2f} (>= 1.7)")
    report("momentum conserved (periodic waves)", worst_drift <= 1e-12,
           f"largest drift of sum(m U) {worst_drift:.1e} of sum(m |U|) (<= 1e-12)")
    report("Lagrangian density rho = m / V_L", worst_rho <= 1e-13,
           f"largest |rho V_L - m| / m {worst_rho:.1e} (<= 1e-13)")

    # First-order switch: runs, converges at about first order
    errs = []
    for n in (50, 100):
        p = wave_case(work, f"first{n}", n, "P", (1, 0, 0), secondOrder=False)
        errs.append(sine_error(p, "P", (1, 0, 0))[0])
    order = math.log2(errs[0] / errs[1])
    report("secondOrder no (first-order switch)", 0.7 <= order <= 1.3,
           f"L1 errors {errs[0]:.3e} {errs[1]:.3e}, order {order:.2f} (0.7 to 1.3)")


def check_stability(work, maxCo):
    """45 degree (2D) and body-diagonal (3D) step P-waves over 2 crossings at
    the default maxCo: no growth (largest |U| in the second crossing not above
    the first), and the 45 degree overshoot."""
    for dim, n, k in [(2, 50, (1, 1, 0)), (3, 24, (1, 1, 1))]:
        p = wave_case(work, f"stable{dim}d", n, "P", k, profile="step",
                      crossings=2.0, dim=dim, maxCo=maxCo)
        h = cs.history(p)
        u = h[:, 9] / (CP * A)
        half = len(u) // 2
        first, second = u[1:half].max(), u[half:].max()
        what = "45 deg" if dim == 2 else "body-diagonal"
        report(f"stable at maxCo {maxCo} in {dim}D ({what} step)",
               second <= first and first < 1.1,
               f"largest |U| / (c_p e) {first:.4f} in the first crossing, "
               f"{second:.4f} in the second (no growth); overshoot "
               f"{100 * (first - 1):.2f} %")


def principal_2d(s):
    """In-plane principal stresses (descending) and the angle (deg) of the
    first, from (xx, xy, yy)."""
    xx, xy, yy = s
    c = 0.5 * (xx + yy)
    r = math.hypot(0.5 * (xx - yy), xy)
    return c + r, c - r, math.degrees(0.5 * math.atan2(2 * xy, xx - yy))


def check_rotation(work):
    """Uniform pre-stress under prescribed rigid rotation, one full turn:
    principal stresses within 0.1 %, direction within 0.1 degree, at every
    step against the exact rotation."""
    s0 = [-100e3, 30e3, 0.0, -50e3, 0.0, -60e3]
    w = 10.0
    check = dict(setup="rotation", omega=[0.0, 0.0, w], centre=[0.5, 0.5, 0.0],
                 turns=1, steps=2000, preStress=s0)
    sides = {s: "fixed" for s in ("xMin", "xMax", "yMin", "yMax")}
    p = cs.run(cs.make_case(os.path.join(work, "rotation"), 8, sides, check,
                            solveMomentum=False))
    h = cs.history(p)
    p1, p2, a0 = principal_2d((s0[0], s0[1], s0[3]))
    scale = max(abs(p1), abs(p2))
    perr, aerr = 0.0, 0.0
    for row in h:
        q1, q2, a = principal_2d((row[2], row[3], row[5]))
        perr = max(perr, abs(q1 - p1) / scale, abs(q2 - p2) / scale,
                   abs(row[7] - s0[5]) / scale)
        da = (a - (a0 + math.degrees(w * row[1])) + 90) % 180 - 90
        aerr = max(aerr, abs(da))
    spread = h[:, 8].max() / scale
    report("rigid rotation, principal stresses", perr <= 1e-3,
           f"largest error {100 * perr:.2e} % over one turn (<= 0.1 %), "
           f"cells differ by {spread:.1e}")
    report("rigid rotation, principal direction", aerr <= 0.1,
           f"largest error {aerr:.2e} deg over one turn (<= 0.1 deg)")


def largest_drop(x):
    """Largest decrease of x after its running maximum (oscillation)."""
    return float(np.max(np.maximum.accumulate(x) - x))


def check_shear(work):
    """Uniform simple shear to gamma 1: von Mises s_xy = Y/sqrt(3) after yield
    within 0.1 %, no oscillation; Drucker-Prager on its cone."""
    sides = {s: "fixed" for s in ("xMin", "xMax", "yMin", "yMax")}
    base = dict(setup="shear", gammaDot=1.0, gammaEnd=1.0, steps=2000,
                centre=[0.5, 0.5, 0.0])

    vm = dict(base, preStress=[0.0] * 6)
    p = cs.run(cs.make_case(os.path.join(work, "shearVM"), 8, sides, vm,
                            material="vonMises", solveMomentum=False))
    h = cs.history(p)
    gam = h[:, 1]
    k = 10e3 / math.sqrt(3)
    after = gam > 2 * k / G
    tau = h[:, 3]
    dev = np.abs(tau[after] - k).max() / k
    drop = largest_drop(tau) / k
    report("simple shear, von Mises", dev <= 1e-3 and drop <= 1e-3,
           f"after yield s_xy within {100 * dev:.2e} % of Y/sqrt(3) (<= 0.1 %), "
           f"largest drop {100 * drop:.2e} % (<= 0.1 %); at gamma 1 s_xx "
           f"{h[-1, 2]:.3f}, s_yy {h[-1, 5]:.3f} Pa")

    p0 = -50e3
    dp = dict(base, preStress=[p0, 0.0, 0.0, p0, 0.0, p0])
    p = cs.run(cs.make_case(os.path.join(work, "shearDP"), 8, sides, dp,
                            material="druckerPrager", solveMomentum=False))
    h = cs.history(p)
    eta, xi, _ = models.dp_coefficients(
        dict(E=cs.E, nu=cs.NU, phi=30.0, c=5e3, psi=0.0, fit="compression"))
    s = h[:, 2:8]       # xx xy xz yy yz zz
    pm = (s[:, 0] + s[:, 3] + s[:, 5]) / 3
    dev = s.copy()
    for i in (0, 3, 5):
        dev[:, i] -= pm
    sqj2 = np.sqrt(0.5 * (dev[:, 0] ** 2 + dev[:, 3] ** 2 + dev[:, 5] ** 2)
                   + dev[:, 1] ** 2 + dev[:, 2] ** 2 + dev[:, 4] ** 2)
    lim = xi * 5e3 - eta * pm
    yielded = h[:, 1] > 2 * lim / G
    err = np.abs(sqj2[yielded] - lim[yielded]).max() / lim[yielded].min()
    pdev = np.abs(pm - p0).max() / abs(p0)
    drop = largest_drop(h[:, 3]) / lim.min()
    report("simple shear, Drucker-Prager", err <= 1e-3 and drop <= 1e-3,
           f"after yield sqrt(J2) within {100 * err:.2e} % of xi c - eta p "
           f"(<= 0.1 %), mean stress change {pdev:.1e}, s_xy largest drop "
           f"{100 * drop:.2e} %")


def check_reflection(work):
    """Step pulses along x reflecting from a boundary: reflected stress and
    velocity, amplitude and sign within 2 %. P-pulse at a fixed wall
    (fixedValue U) and a free surface (zeroGradient U); P- and S-pulse at a
    symmetry plane (no normal velocity: like a wall for P; no shear
    traction: like a free surface for S)."""
    n = 200
    for kind, mode in [("fixed", "P"), ("free", "P"), ("symmetry", "P"),
                       ("symmetry", "S")]:
        sides = {"xMin": "free", "xMax": kind, "yMin": "empty", "yMax": "empty"}
        p = wave_case(work, f"reflect_{kind}{mode}", n, mode, (1, 0, 0),
                      profile="pulse", crossings=0.6, sides=sides, ny=1,
                      pulse=(0.5, 0.7))
        c = cs.cells(p)
        x = c[:, 0]
        mid = (x > 0.75) & (x < 0.85)
        if mode == "P":
            # incident: s_xx = M e, U_x = -c_p e
            s = c[mid, 6].mean() / (M * A)
            u = c[mid, 3].mean() / (-CP * A)
            behind = np.abs(c[x < 0.6, 6]).max() / (M * A)
        else:
            # incident: s_xy = G e, U_y = -c_s e. Run time is set from c_s,
            # so the reflected pulse is at the same place.
            s = c[mid, 7].mean() / (G * A)
            u = c[mid, 4].mean() / (-CS * A)
            behind = np.abs(c[x < 0.6, 7]).max() / (G * A)
        # From a wall: same stress, opposite velocity; from a free surface:
        # opposite stress, same velocity
        wall = kind == "fixed" or (kind == "symmetry" and mode == "P")
        want = (1.0, -1.0) if wall else (-1.0, 1.0)
        err = max(abs(s - want[0]), abs(u - want[1]))
        name = {"fixed": "fixed wall", "free": "free surface",
                "symmetry": "symmetry plane"}[kind]
        report(f"{mode}-pulse reflection from a {name}", err <= 0.02,
               f"reflected stress {s:+.4f}, velocity {u:+.4f} of the incident "
               f"(want {want[0]:+.0f}, {want[1]:+.0f}, within 2 %); behind it "
               f"|stress| <= {behind:.1e}")


def check_parallel(work):
    """45 degree wave on 2 and 4 processors matches serial to 1e-12."""
    def load(path, nProcs):
        c = cs.cells(path, nProcs)
        order = np.lexsort(np.round(c[:, 2::-1].T, 9))
        return c[order]

    ref = None
    for nProcs, split in [(1, None), (2, (2, 1, 1)), (4, (2, 2, 1))]:
        p = wave_case(work, f"par{nProcs}", 50, "P", (1, 1, 0), profile="step",
                      crossings=1.0, nProcs=nProcs, split=split)
        c = load(p, nProcs)
        if ref is None:
            ref = c
            continue
        du = np.abs(c[:, 3:6] - ref[:, 3:6]).max() / (CP * A)
        ds = np.abs(c[:, 6:12] - ref[:, 6:12]).max() / (M * A)
        drift = momentum_drift(p)
        report(f"{nProcs} processors match serial",
               c.shape == ref.shape and max(du, ds) <= 1e-12 and drift <= 1e-12,
               f"largest difference U {du:.1e}, sigma {ds:.1e} (<= 1e-12); "
               f"momentum drift {drift:.1e}")


def main():
    work = sys.argv[1]
    maxCo = float(sys.argv[2])
    check_waves(work)
    check_stability(work, maxCo)
    check_rotation(work)
    check_shear(work)
    check_reflection(work)
    check_parallel(work)
    t = cs.log_value(os.path.join(work, "waveP0_100"), "Run time per cell per step")
    print(f"run time per cell per step (P-wave, 100 x 100, serial): {t:.2e} s")
    if failures:
        print(f"{len(failures)} check(s) failed: {', '.join(failures)}")
        sys.exit(1)
    print("all checks passed")


if __name__ == "__main__":
    main()
