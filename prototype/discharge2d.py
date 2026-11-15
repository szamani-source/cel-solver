#!/usr/bin/env python3
"""WP 3.3p: target T0 in the 2D prototype, settle then discharge.

A rectangular tank (docs/spec/target-problems.md, T0) is filled with
cohesionless Mohr-Coulomb material, settles under gravity with the outlet
closed (WP 3.1p controls: gravity ramp 0.1 s, damping 250 1/s, decision 0020),
then the outlet, a centred segment of the bottom wall, becomes an outflow face
(rule 10 of ``void2d.py``: zero traction, empty ghosts, material swept across
it leaves the domain and is counted). Damping is off during the discharge.

The free surface, the void above the material and the face shear cap come
from ``void2d.VoidBox`` (WP 3.2p, 3.2q; decisions 0021, 0023).

``python3 discharge2d.py --check`` runs the outlet checks
(tests/test_proto_silo_discharge.sh; 10 mm cells); ``--rate-check`` the rate
law checks (WP 3.3q, decision 0027; 5 mm cells, four processes). Both at
E = 0.1 MPa.
``--run`` runs one case and saves its history (``--save`` .npz) for
``tools/plot_discharge.py``. SI units, thickness 1 m.
"""

import argparse
import math
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import silo2d as sd  # noqa: E402
import void2d as vd  # noqa: E402

W_TANK = 0.40
H_FILL = 0.60
D_OUT = 0.06
DELTA_W = 20.0
T_RAMP = 0.1
ALPHA = 250.0
T_SETTLE = 0.3           # ramp 0.1 s plus 0.2 s damped (decision 0020)
HEADROOM = 0.02          # void above the fill (m); the material never rises
WINDOW = (0.95, 0.70)    # steady-rate window, fractions of the initial mass
T_AVG = 0.3              # flow pattern: velocities averaged over this much more discharge,
                         # at least one period of the bed's first elastic mode (4 H / c_p)


def make_tank(cell=0.005, W=W_TANK, H=H_FILL, D=D_OUT, E=None, delta_w=DELTA_W,
              low_mach=False):
    """Tank of width W with material to height H, void above, frictional
    walls and floor (delta_w), outlet of width D centred in the floor
    (closed until opened)."""
    nx = int(round(W / cell))
    ny = int(math.ceil((H + HEADROOM) / cell))
    params = dict(sd.T0_MC)
    if E is not None:
        params["E"] = E
    walls = {"left": "friction", "right": "friction", "bottom": "friction", "top": "smooth"}
    box = vd.make_box(nx, ny, nx * cell, ny * cell, params=params, walls=walls,
                      delta_w=delta_w)
    box.fill(vd.rect(-1.0, 2.0 * W, -1.0, H, cell, cell))     # stress free
    box.t_ramp = T_RAMP
    box.low_mach = False          # settling keeps the full jump terms (rule 11 off)
    box.low_mach_open = low_mach  # switched on with the outlet
    box.damping = ALPHA
    n_out = int(round(D / cell))
    if abs(n_out * cell - D) > 1e-9 or (nx - n_out) % 2:
        raise ValueError(f"outlet {D} m is not a centred whole number of {cell} m cells")
    j0 = (nx - n_out) // 2
    box.outlet_cols = np.zeros(nx, dtype=bool)
    box.outlet_cols[j0:j0 + n_out] = True
    box.D = D
    box.H_fill = H
    return box


def settle(box, ke_tol=1e-4, t_settle=T_SETTLE, cfl=0.5):
    """Stage 1: gravity ramp plus damping for t_settle. The walls sit at the
    friction limit and the column creeps slowly (WP 3.1p), so a fixed time is
    used and the kinetic energy is checked against ke_tol times the strain
    energy. Returns (settled, wall seconds, KE / strain energy)."""
    t0 = time.perf_counter()
    while box.t < t_settle - 1e-12:
        box.step(min(box.stable_dt(cfl), t_settle - box.t))
    ratio = box.kinetic_energy() / box.strain_energy()
    return ratio <= ke_tol, time.perf_counter() - t0, ratio


def discharge(box, t_end=None, m_stop=0.5, cfl=0.5, every=None, frames=None,
              frame_dt=0.02, log=None):
    """Stage 2: open the outlet at the current time and run without damping
    until t_end (s after opening) or until the mass in the tank is m_stop
    times the mass at opening. Records (t - t_open, mass in tank, discharged)
    every `every` steps (default: every 2 ms of time) and, if `frames` is a
    list, appends snapshots every frame_dt."""
    box.damping = 0.0
    box.low_mach = getattr(box, "low_mach_open", False)
    box.outlet = box.outlet_cols
    box.t_open = box.t
    m0 = box.total_mass()
    box.m_open = m0
    d0 = box.discharged
    hist = [(0.0, m0, 0.0)]
    t_rec = 0.002
    next_rec = t_rec
    next_frame = 0.0
    t0 = time.perf_counter()
    while True:
        tau = box.t - box.t_open
        if frames is not None and tau >= next_frame - 1e-12:
            frames.append(snapshot(box))
            next_frame += frame_dt
        if (t_end is not None and tau >= t_end) or box.total_mass() <= m_stop * m0:
            break
        box.step(box.stable_dt(cfl))
        tau = box.t - box.t_open
        if tau >= next_rec:
            hist.append((tau, box.total_mass(), box.discharged - d0))
            next_rec += t_rec
            if log:
                log(box, hist[-1])
    box.wall_seconds = time.perf_counter() - t0
    return np.array(hist)


def snapshot(box):
    """Fields for the figures and the animation."""
    return dict(t=box.t - getattr(box, "t_open", box.t), al=box.al.copy(),
                vx=box.vx.copy(), vy=box.vy.copy(), epq=box.epq.copy(),
                p=-(box.s[..., 0] + box.s[..., 1] + box.s[..., 2]) / 3.0 * box.al,
                mass=box.total_mass(), out=box.discharged)


# ----------------------------------------------------------------------------
# Measurements
# ----------------------------------------------------------------------------

def steady_rate(hist, m_open, window=WINDOW):
    """Steady discharge rate (kg/s/m): least-squares slope of the tank mass
    against time while it falls from window[0] to window[1] times m_open.
    Also returns the largest deviation of 0.1 s sub-window rates from it."""
    t, m = hist[:, 0], hist[:, 1]
    sel = (m <= window[0] * m_open) & (m >= window[1] * m_open)
    if sel.sum() < 5:
        return float("nan"), float("nan"), (float("nan"), float("nan"))
    ts, ms = t[sel], m[sel]
    rate = -np.polyfit(ts, ms, 1)[0]
    subs = []
    span = 0.1
    a = ts[0]
    while a + span <= ts[-1] + 1e-9:
        k = (ts >= a) & (ts <= a + span)
        if k.sum() >= 3:
            subs.append(-np.polyfit(ts[k], ms[k], 1)[0])
        a += span
    var = max(abs(r / rate - 1.0) for r in subs) if subs else float("nan")
    return rate, var, (ts[0], ts[-1])


def beverloo_C(rate, D, rho=sd.RHO, g=sd.G_ACC):
    return rate / (rho * math.sqrt(g) * D ** 1.5)


def flow_pattern(box, corner=0.05, vx=None, vy=None):
    """Largest speed in the two bottom corner squares (side `corner`) over the
    mean downward speed of the material in the outlet cells, and the crater
    depth (column height at the walls minus at the centre). vx, vy: velocity
    fields to use instead of the box's (time averages, averaged_pattern)."""
    vx = box.vx if vx is None else vx
    vy = box.vy if vy is None else vy
    sp = np.hypot(vx, vy)
    X, Y = np.meshgrid(box.xc, box.yc)
    has = box.al > 0.5
    cz = has & (Y < corner) & ((X < corner) | (X > box.lx - corner))
    v_out = -(vy[0, box.outlet_cols] * box.rho[0, box.outlet_cols]).sum() / max(
        box.rho[0, box.outlet_cols].sum(), 1e-300)
    hc = box.column_height()
    n = len(hc)
    crater = 0.5 * (hc[:3].mean() + hc[-3:].mean()) - hc[n // 2 - 1:n // 2 + 1].mean()
    return dict(corner_ratio=(sp[cz].max() if cz.any() else 0.0) / max(v_out, 1e-300),
                v_out=v_out, crater=crater)


def averaged_pattern(box, t_avg=T_AVG, cfl=0.5):
    """flow_pattern with the velocities averaged over a further t_avg of
    discharge (WP 3.3q, decision 0027). A soft bed (E = 0.1 MPa) rings
    elastically after the outlet opens: the instantaneous corner speed is 2 to
    6 percent of the outlet speed while the corner hardly moves on average
    (0.6 percent over 0.4 s). The stagnant-corner check is about flow, so it
    uses the mean velocity over at least one period of the bed's first
    elastic mode (4 H / c_p = 0.28 s for H = 0.6 m at E = 0.1 MPa; 0.1 s
    still left 1.5 percent). The instantaneous ratio is kept as
    "corner_ratio_now"."""
    now = flow_pattern(box)["corner_ratio"]
    ax = np.zeros_like(box.vx)
    ay = np.zeros_like(box.vy)
    T = 0.0
    while T < t_avg - 1e-12:
        dt = min(box.stable_dt(cfl), t_avg - T)
        box.step(dt)
        ax += dt * box.vx
        ay += dt * box.vy
        T += dt
    pat = flow_pattern(box, vx=ax / T, vy=ay / T)
    pat["corner_ratio_now"] = now
    return pat


def run_case(cell=0.005, W=W_TANK, H=H_FILL, D=D_OUT, E=None, window=WINDOW, frames=None,
             frame_dt=0.02, m_stop=None, t_end=None, log=None, quiet=False, low_mach=False,
             average=True):
    """Settle, open, discharge to the end of the rate window. Returns a dict."""
    box = make_tank(cell, W, H, D, E, low_mach=low_mach)
    m_fill = box.total_mass()
    ok, sec_settle, ke_ratio = settle(box)
    t_settle = box.t
    m_open = box.total_mass()
    leaked = box.discharged
    if not quiet:
        print(f"    settled {ok} (KE / strain energy {ke_ratio:.1e}) at t = {t_settle:.3f} s ({box.nstep} steps, "
              f"{sec_settle:.0f} s); mass through the closed outlet {leaked:.1e} kg/m",
              flush=True)
    stop = m_stop if m_stop is not None else window[1] - 0.02
    hist = discharge(box, t_end=t_end, m_stop=stop, frames=frames, frame_dt=frame_dt,
                     log=log)
    t_stop = box.t - box.t_open
    pat = averaged_pattern(box) if average else flow_pattern(box)
    pat["t"] = t_stop
    rate, var, win = steady_rate(hist, m_open, window)
    cons = np.abs(hist[:, 1] + hist[:, 2] - m_open).max() / m_open
    return dict(box=box, hist=hist, settled=ok, ke_ratio=ke_ratio, t_settle=t_settle, leaked=leaked,
                m_fill=m_fill, m_open=m_open, rate=rate, var=var, win=win,
                C=beverloo_C(rate, D), cons=cons, pattern=pat,
                seconds=sec_settle + box.wall_seconds, steps=box.nstep, cell=cell,
                W=W, H=H, D=D, E=box.params["E"])


def fit_exponent(Ds, rates):
    """Free power-law fit rate = a D^n, and C with n fixed at 1.5."""
    n, lna = np.polyfit(np.log(Ds), np.log(rates), 1)
    Cs = [beverloo_C(r, D) for r, D in zip(rates, Ds)]
    # least squares C for rate = C rho sqrt(g) D^1.5
    x = np.array([sd.RHO * math.sqrt(sd.G_ACC) * D ** 1.5 for D in Ds])
    C = float((x * np.array(rates)).sum() / (x * x).sum())
    return n, C, Cs


# ----------------------------------------------------------------------------
# CI checks (coarse)
# ----------------------------------------------------------------------------

CHECK_CELL = 0.01
RATE_CELL = 0.005        # rate-law checks (WP 3.3q, decision 0027): D = 0.04 m is 8 cells
CHECK_E = 1.0e5
CHECK_WINDOW = (0.97, 0.88)
CHECK_C_MIN = 0.8        # guard for rule 10b (decision 0027): 0.98 with it, about 0.3 without


def _check_case(args):
    """One check case (run in a worker process); the box is not returned."""
    cell, H, D = args
    t0 = time.perf_counter()
    r = run_case(cell=cell, H=H, D=D, E=CHECK_E, window=CHECK_WINDOW, quiet=True)
    r.pop("box")
    r.pop("hist")
    r["wall"] = time.perf_counter() - t0
    return r


def check(full=False, procs=None):
    """WP 3.3p and 3.3q checks, E = 0.1 MPa (the stiffest E at which the
    checks run in minutes; the rate still depends on E, decision 0027), rate
    measured from 97 to 88 percent of the mass. Cases run in parallel
    processes (procs, default up to 4).

    Default (CI, tests/test_proto_silo_discharge.sh), 10 mm cells: the outlet
    mechanics, D = 0.06 m at H = 0.6 and 0.9 m: nothing leaves before opening,
    mass in tank plus discharged constant to 1e-8, head independence within
    10 %, and C at least CHECK_C_MIN (the outlet is not held shut, rule 10b).
    full=True (tests/test_proto_silo_rate.sh, WP 3.3q), 5 mm cells so that
    the narrowest outlet is 8 cells wide (at 10 mm, 4 cells, it gives C 0.75
    against 0.98 for the others): also D = 0.04 and 0.08 m, steady rate within
    10 %, exponent 1.5 +- 0.2, C in 1.0 to 1.9, stagnant corners (mean speed
    over 0.3 s) below 1 % of the outlet speed, a crater."""
    import multiprocessing as mp
    ok = True
    t_all = time.perf_counter()
    P = lambda good: "PASS" if good else "FAIL"            # noqa: E731
    cell = RATE_CELL if full else CHECK_CELL
    print(f"T0 tank W = {W_TANK} m, frictional walls and floor ({DELTA_W:g} deg), "
          f"Mohr-Coulomb c = 0, phi 35, E = {CHECK_E / 1e6:g} MPa, {1000 * cell:g} mm cells")
    cases = [(0.6, 0.06), (0.9, 0.06)]
    if full:
        cases = [(0.9, 0.06), (0.6, 0.04), (0.6, 0.06), (0.6, 0.08)]   # longest first
    n = procs or min(4, os.cpu_count() or 1, len(cases))
    print(f"{len(cases)} cases on {n} processes", flush=True)
    with mp.get_context("fork").Pool(n) as pool:
        out = pool.map(_check_case, [(cell, H, D) for H, D in cases], chunksize=1)
    res = dict(zip(cases, out))
    for H, D in sorted(cases, key=lambda c: (c[0], c[1])):
        r = res[(H, D)]
        print(f"\n[H = {H} m, D = {D} m]")
        print(f"    settled at t = {r['t_settle']:.3f} s (KE / strain energy {r['ke_ratio']:.1e})")
        g_leak = r["leaked"] == 0.0
        g_cons = r["cons"] <= 1e-8
        g_flow = r["rate"] > 0.0 and r["pattern"]["v_out"] > 0.0
        g_var = r["var"] <= 0.10
        good = g_leak and g_cons and g_flow and (g_var or not full)
        ok &= good
        print(f"  rate {r['rate']:.1f} kg/s/m over t = {r['win'][0]:.3f} to "
              f"{r['win'][1]:.3f} s after opening (C = {r['C']:.3f}); 0.1 s windows within "
              f"{100 * r['var']:.1f} %{' (<= 10 %)' if full else ''}; mass in tank + "
              f"discharged constant to {r['cons']:.0e} (<= 1e-8); mass through the closed "
              f"outlet {r['leaked']:.0e} (0); {r['steps']} steps, {r['wall']:.0f} s: "
              f"{P(good)}")
        pat = r["pattern"]
        print(f"  at the end (t = {pat['t']:.2f} s, velocities averaged over a further "
              f"{T_AVG:g} s): outlet speed {pat['v_out']:.2f} m/s, corner speed / outlet speed "
              f"{pat['corner_ratio']:.1e} (instantaneous {pat['corner_ratio_now']:.1e}), crater "
              f"depth {1000 * pat['crater']:.0f} mm")

    r6 = res[(0.6, 0.06)]
    g = r6["C"] >= CHECK_C_MIN
    ok &= g
    print(f"\n[outlet open] D = 0.06 m, H = 0.6 m: C = {r6['C']:.3f} (>= {CHECK_C_MIN}; 0.3 "
          f"before rule 10b, decision 0027): {P(g)}")
    rh = res[(0.9, 0.06)]["rate"] / res[(0.6, 0.06)]["rate"]
    g = abs(rh - 1.0) <= 0.10
    ok &= g
    print(f"[head] rate at H = 0.9 m / at H = 0.6 m = {rh:.3f} (within 10 %): {P(g)}")
    if full:
        Ds = [0.04, 0.06, 0.08]
        rates = [res[(0.6, D)]["rate"] for D in Ds]
        n, C, Cs = fit_exponent(Ds, rates)
        g = abs(n - 1.5) <= 0.2 and 1.0 <= C <= 1.9
        ok &= g
        print(f"[rate law] exponent {n:.3f} (1.3 to 1.7), C = {C:.3f} with exponent 1.5 "
              f"(1.0 to 1.9; reference 1.45), per outlet "
              f"{', '.join(f'{c:.3f}' for c in Cs)}: {P(g)}")
        pat = res[(0.6, 0.06)]["pattern"]
        g = pat["corner_ratio"] <= 0.01 and pat["crater"] > 0.0
        ok &= g
        print(f"[pattern] D = 0.06 m: corner speed {pat['corner_ratio']:.1e} of the outlet "
              f"speed (<= 1e-2), crater {1000 * pat['crater']:.0f} mm (> 0): {P(g)}")
    print(f"\nWP 3.3{'q rate' if full else 'p outlet'} checks: {'ALL PASS' if ok else 'FAILED'} "
          f"({time.perf_counter() - t_all:.0f} s)")
    return ok


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--check", action="store_true", help="run the CI checks (outlet)")
    p.add_argument("--rate-check", action="store_true",
                   help="run the rate-law checks (tests/test_proto_silo_rate.sh, WP 3.3q)")
    p.add_argument("--run", action="store_true", help="run one case")
    p.add_argument("--cell", type=float, default=0.005)
    p.add_argument("--H", type=float, default=H_FILL)
    p.add_argument("--D", type=float, default=D_OUT)
    p.add_argument("--E", type=float, default=None)
    p.add_argument("--low-mach", type=float, default=None,
                   help="rule 11 on after opening, with this floor (study; noisy)")
    p.add_argument("--frames", action="store_true", help="store snapshots for the animation")
    p.add_argument("--frame-dt", type=float, default=0.02)
    p.add_argument("--m-stop", type=float, default=None,
                   help="stop at this fraction of the mass at opening")
    p.add_argument("--t-end", type=float, default=None, help="stop this long after opening")
    p.add_argument("--save", default=None, help=".npz file for the history and frames")
    args = p.parse_args(argv)
    if args.check or args.rate_check:
        return 0 if check(full=args.rate_check) else 1
    if args.run:
        frames = [] if args.frames else None
        t0 = time.perf_counter()

        def log(b, rec):
            if int(round(rec[0] / 0.002)) % 50 == 0:
                print(f"    t = {rec[0]:.2f} s  tank {rec[1]:.1f} kg/m  out {rec[2]:.1f} kg/m  "
                      f"({time.perf_counter() - t0:.0f} s)", flush=True)
        r = run_case(cell=args.cell, H=args.H, D=args.D, E=args.E, frames=frames,
                     frame_dt=args.frame_dt, m_stop=args.m_stop, log=log, t_end=args.t_end,
                     low_mach=False if args.low_mach is None else (args.low_mach or True))
        print(f"rate {r['rate']:.2f} kg/s/m, C = {r['C']:.4f}, variation {100 * r['var']:.1f} %, "
              f"window {r['win'][0]:.3f} to {r['win'][1]:.3f} s, conservation {r['cons']:.0e}, "
              f"leak {r['leaked']:.0e}, pattern {r['pattern']}, {r['steps']} steps, "
              f"{r['seconds']:.0f} s")
        if args.save:
            b = r["box"]
            extra = {}
            if frames:
                for k in ("al", "vx", "vy", "epq", "p"):
                    extra["f_" + k] = np.array([f[k] for f in frames])
                extra["f_t"] = np.array([f["t"] for f in frames])
            np.savez_compressed(
                args.save, hist=r["hist"], m_open=r["m_open"], rate=r["rate"], C=r["C"],
                var=r["var"], win=np.array(r["win"]), cell=r["cell"], W=r["W"], H=r["H"],
                D=r["D"], E=r["E"], xc=b.xc, yc=b.yc, lx=b.lx, ly=b.ly,
                outlet=b.outlet_cols, seconds=r["seconds"], t_settle=r["t_settle"],
                cons=r["cons"], **extra)
        return 0
    p.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
