"""WP F.3 done-when: gas-filled column loaded at the top (decisions 0125, 0128).

Two-way coupling: the solid's sigma is the effective stress, the gas
pressure acts on the material, and the grains carry and squeeze the gas.
Every run uses the ergunWenYu drag law. Column of tests/gas/
consolidation_cases.py: fine quartz powder (d 0.1 mm, gas fraction 0.40,
E 1 MPa, nu 0.3), 1 m tall, one cell wide, smooth sealed walls and base,
drained top loaded with 2 kPa at t = 0, velocity damping 50 1/s.

A  Consolidation: 50 cells, to T = 1 (1.15 s). Degree of consolidation
   U = 1 - mean(p)/p_i against Terzaghi at T = 0.05, 0.1, 0.2, 0.5 and 1
   within 3 %, with c = k/(mu (e/p0 + 1/M)) and p_i = q (1/M)/(e/p0 + 1/M).
   (Earlier than T = 0.05, 0.36 s, the elastic waves of the sudden load,
   left out by the quasi-static theory, are still being damped.) The
   pressure profile is reported against the series. 0.1 mm grains at the
   v1 step: the steps are those of the same column without gas (number
   and smallest step within 1 %), no blow-up.
R  At rest: the column under gravity and the load, started in its final
   equilibrium (effective stress -(q + rho g depth), K0 sideways), run for
   0.5 s: the pore pressure stays at p0 to 1e-6 relative.
P  A to T = 0.2 on 2 and 4 processors (cuts across the column): p and U
   match serial to 1e-10 (of p_i and of the largest speed); gas mass less
   the boundary inflow is constant to 1e-12 on every count.
C  Cost per step with gas on against gas off: the 1D column (A) and a
   20 x 100 cell column to T = 0.05 (report only; F.9 sets the 1.5 limit).

    python3 consolidation_checks.py WORKDIR
"""

import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import consolidation_cases as cc  # noqa: E402
import pore_cases as pc  # noqa: E402

TS = (0.05, 0.1, 0.2, 0.5, 1.0)
NY = 50
res = []


def report(name, good, detail):
    res.append(bool(good))
    print(f"  {name}: {detail}: {'PASS' if good else 'FAIL'}", flush=True)


def at_T(path):
    """Write times keyed by the nearest T in TS."""
    out = {}
    for t in pc.times(path):
        T = float(t) * cc.C_V / cc.H ** 2
        for x in TS:
            if abs(T - x) < 1e-6:
                out[x] = t
    return out


def drift(path):
    m, ni, _, _ = pc.mass_log(path)
    return float(np.abs((m - ni) - (m[0] - ni[0])).max() / m[0])


def main(work):
    print(f"Column: M {cc.M_MOD:.6g} Pa, c {cc.C_V:.6g} m^2/s, p_i {cc.P_I:.6g} Pa "
          f"({cc.P_I / cc.Q:.4f} of the load), T = 1 at {cc.t_of(1):.4g} s; damping "
          f"delay estimate rho D c / M = {cc.RHO * cc.DAMPING * cc.C_V / cc.M_MOD:.2%}")

    # A
    print("A: consolidation against Terzaghi (50 cells)")
    a = f"{work}/A"
    cc.make(a, ny=NY, end=cc.t_of(1.0), write_interval=cc.t_of(0.05))
    wall_a = pc.run(a)
    ts = at_T(a)
    for T in TS:
        U = cc.degree(a, ts[T])
        Ut = cc.theory(T)
        prof = pc.field(a, ts[T], "p") / cc.P_I
        perr = float(np.abs(prof - cc.profile_theory(NY, T)).max())
        report(f"T = {T}", abs(U - Ut) <= 0.03 * Ut,
               f"U {U:.4f} against {Ut:.4f} ({(U - Ut) / Ut:+.2%}); profile within "
               f"{perr:.3f} p_i")
    off = f"{work}/A_off"
    cc.make(off, ny=NY, end=cc.t_of(1.0), write_interval=cc.t_of(1.0), gas=False)
    wall_off = pc.run(off)
    dt_on, dt_off = cc.deltaT(a), cc.deltaT(off)
    # The steps follow the density (sqrt(M/rho)): with gas the column
    # compresses later, so they differ slightly; the gas sets no limit
    nr = len(dt_on) / len(dt_off) - 1
    mr = dt_on.min() / dt_off.min() - 1
    same = abs(nr) < 0.01 and abs(mr) < 0.01
    pmax = max(float(np.abs(pc.field(a, t, "p")).max()) for t in pc.times(a))
    report("0.1 mm grains at the v1 step", same and np.isfinite(pmax) and pmax < 1.5 * cc.P_I,
           f"{len(dt_on)} steps with gas, {len(dt_off)} without ({nr:+.1e}), smallest step "
           f"{dt_on.min():.4e} s against {dt_off.min():.4e} s; largest |p| {pmax:.1f} Pa")
    report("gas mass balance (A)", drift(a) <= 1e-12, f"{drift(a):.1e}")

    # R
    print("R: column at rest under gravity and the load")
    r = f"{work}/R"
    cc.make(r, ny=NY, end=0.5, write_interval=0.1, at_rest=True, gravity=9.81)
    pc.run(r)
    worst = max(float(np.abs(pc.field(r, t, "p")).max()) for t in pc.times(r)) / pc.P0
    umax = max(float(np.abs(pc.field(r, t, "U", NY)).max()) for t in pc.times(r))
    report("pore pressure at rest", worst <= 1e-6,
           f"largest |p - p0|/p0 {worst:.1e} over 0.5 s (speeds up to {umax:.1e} m/s)")

    # P
    print("P: 2 and 4 processors against serial (to T = 0.2)")
    runs = {}
    secs = {}
    for n, dec in ((1, (1, 1, 1)), (2, (1, 2, 1)), (4, (1, 4, 1))):
        w = f"{work}/P{n}"
        cc.make(w, ny=NY, end=cc.t_of(0.2), write_interval=cc.t_of(0.05), decomp=dec)
        secs[n] = pc.run(w, n)
        runs[n] = w
        report(f"gas mass balance, {n} processor(s)", drift(w) <= 1e-12, f"{drift(w):.1e}")
    for n in (2, 4):
        dp = du = 0.0
        for t in pc.times(runs[1]):
            dp = max(dp, float(np.abs(pc.field(runs[1], t, "p") - pc.field(runs[n], t, "p")).max()))
            u1 = pc.field(runs[1], t, "U", NY)
            un = pc.field(runs[n], t, "U", NY)
            du = max(du, float(np.abs(u1 - un).max() / max(np.abs(u1).max(), 1e-300)))
        report(f"{n} processors against serial", dp / cc.P_I <= 1e-10 and du <= 1e-10,
               f"p within {dp / cc.P_I:.1e} of p_i, U within {du:.1e}")
    print(f"  wall times: 1 proc {secs[1]:.2f} s, 2 procs {secs[2]:.2f} s, 4 procs {secs[4]:.2f} s "
          f"({NY} cells: too small to scale)")

    # C
    print("C: cost per step, gas on against gas off (report)")
    print(f"  1D column (A, {NY} cells, to T = 1): {wall_a:.2f} s against {wall_off:.2f} s: "
          f"{wall_a / wall_off:.2f}x")
    w2 = {}
    for gas in (False, True):
        w = f"{work}/C{int(gas)}"
        cc.make(w, nx=20, ny=100, end=cc.t_of(0.05), write_interval=cc.t_of(0.05), gas=gas)
        w2[gas] = (pc.run(w), len(cc.deltaT(w)))
    _, _, solves, _ = pc.mass_log(f"{work}/C1")
    print(f"  20 x 100 cells to T = 0.05 ({w2[True][1]} steps): {w2[True][0]:.2f} s against "
          f"{w2[False][0]:.2f} s: {w2[True][0] / w2[False][0]:.2f}x, "
          f"{solves.mean():.2f} pressure solves per step")

    print(f"{sum(res)}/{len(res)} checks passed")
    return all(res)


if __name__ == "__main__":
    work = sys.argv[1]
    os.makedirs(work, exist_ok=True)
    for d in os.listdir(work):
        shutil.rmtree(os.path.join(work, d), ignore_errors=True)
    sys.exit(0 if main(work) else 1)
