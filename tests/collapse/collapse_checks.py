"""WP 3.2b done-when: granular column collapse in celFoam's cycle with the
face shear cap (rule 8, decision 0023) and top layers spread sideways
(rule 9, decision 0024), neither rule 12 nor 13 (decision 0029); C++ in
decision 0032. Same set-up as the prototype (prototype/void2d.py
column_collapse, tests/test_proto_collapse.sh):

T0 Mohr-Coulomb (rho 1800, nu 0.3, phi 35, c 0, psi 0, no tension), 10 mm
cells (10 across the short side), column L_i x H_i against the smooth left
wall, floor friction 35 deg, released from K0 rest under full gravity, no
damping. a = H_i/L_i = 0.5 at E = 0.1, 1 and 10 MPa, a = 2 at 1 MPa.

Measured as in the prototype: run-out (L_f - L_i)/L_i with L_f the right
edge of the last column holding at least half a cell of mass, final height
H_f (largest column height, sum of alpha dy), deposit slope as a line fitted
where the surface is between 10 and 90 % of H_f (decision 0025), and whether
the collapse stops (kinetic energy below 1e-3 of its peak before t_max, and
the run-out at 1 s against the end).

Checks (decision 0019 bands, decision 0032): every run keeps mass to 1e-10
and leaves no stress in empty cells; at E = 1 MPa (the prototype's CI
setting) a = 0.5 runs out 0.30 to 0.90 with the top kept within a cell and
the fitted slope at most 37 deg, a = 2 runs out 1.2 to 3.6 with the slope at
most 37 deg; at 0.1 and 10 MPa the run-out is recorded against the same band
and the prototype (0.55 and 0.25, decision 0029), not checked.

Run by tests/test_collapse.sh after celVoidCheck is built:

    python3 collapse_checks.py WORKDIR [a:E_MPa ...]
"""

import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "void"))
import void_checks as vc  # noqa: E402  (make_case, run, CELL)

RHO, NU, PHI, DELTA_B = 1800.0, 0.3, 35.0, 35.0
CELL = vc.CELL
T_MAX = 2.0
K0 = NU / (1.0 - NU)
# prototype run-outs, neither rule 12 nor 13 (decision 0029 table; 1 MPa with
# rule 9, decisions 0024 and 0025)
PROTO = {(0.5, 0.1): 0.55, (0.5, 1.0): 0.40, (0.5, 10.0): 0.25, (2.0, 1.0): 1.8}
BAND = {0.5: (0.30, 0.90), 2.0: (1.2, 3.6)}
# WP 3.6e (decision 0050): check every case, not only 1 MPa
CHECK_ALL = os.environ.get("CEL_COLLAPSE_CHECK_ALL", "") == "1"
# history columns (voidCheck)
STEP, T, MASS, SEMPTY, KE, NCAP = 0, 1, 2, 8, 12, 13


def geometry(a):
    if a <= 1.0:
        H_i = 0.1
        L_i = H_i / a
    else:
        L_i = 0.1
        H_i = a * L_i
    run = 1.2 * a * 1.6                       # room for 1.5 times the fit
    nx = int(math.ceil((L_i * (1.0 + run) + 4 * CELL) / CELL - 1e-9))
    ny = int(math.ceil((H_i + 3 * CELL) / CELL - 1e-9))
    return L_i, H_i, nx, ny


def measure(cells, nx, ny, L_i):
    """Run-out, final height and fitted slope from a cells file."""
    c = np.loadtxt(cells)
    i = np.floor(c[:, 0] / CELL + 1e-9).astype(int)
    j = np.floor(c[:, 1] / CELL + 1e-9).astype(int)
    al = np.zeros((ny, nx))
    rho = np.zeros((ny, nx))
    al[j, i] = c[:, 2]
    rho[j, i] = c[:, 3]
    hc = al.sum(axis=0) * CELL
    mcol = rho.sum(axis=0) / RHO               # cells' worth of mass
    front = np.where(mcol >= 0.5)[0]
    L_f = (front.max() + 1) * CELL if front.size else 0.0
    H_f = hc.max()
    band = (hc >= 0.1 * H_f) & (hc <= 0.9 * H_f)
    xc = (np.arange(nx) + 0.5) * CELL
    slope = (math.degrees(math.atan(-np.polyfit(xc[band], hc[band], 1)[0]))
             if band.sum() >= 3 else 0.0)
    return dict(runout=(L_f - L_i) / L_i, H_f=H_f, slope=slope, hc=hc)


def setup(work, a, E_MPa, t_max=T_MAX, name=None, nProcs=1):
    L_i, H_i, nx, ny = geometry(a)
    path = f"{work}/{name or f'collapse_a{a:g}_E{E_MPa:g}'}"
    E = E_MPa * 1.0e6
    material = (f"type mohrCoulomb;\nE {E!r};\nnu {NU!r};\nphi {PHI!r};\nc 0;\npsi 0;\n"
                "sigmaT 0;\nshearCap yes;\nspreadLayers yes;\n")
    check = dict(fill=[0.0, 0.0, L_i, H_i], U0=[0.0, 0.0, 0.0], initK0=K0,
                 tEnd=t_max, maxCo=0.4, writeTimes=[0.5, 1.0])
    vc.make_case(path, nx, ny, material, check, nProcs=nProcs)
    # walls: smooth sides and top, frictional floor (prototype column_collapse)
    mp = f"{path}/constant/materialProperties"
    txt = open(mp).read()
    txt = txt.replace("    bottom { friction 0; }", f"    bottom {{ friction {DELTA_B!r}; }}")
    txt = txt.replace(f"rho {vc.RHO};", f"rho {RHO};")
    open(mp, "w").write(txt)
    return path


def collapse(work, a, E_MPa, t_max=T_MAX):
    L_i, H_i, nx, ny = geometry(a)
    path = setup(work, a, E_MPa, t_max)
    sec = vc.run(path)
    h = np.loadtxt(f"{path}/history.dat")
    ke = h[:, KE]
    peak = ke.max()
    i_peak = int(ke.argmax())
    low = np.where((h[:, T] > 0.05) & (np.arange(len(ke)) > i_peak)
                   & (ke <= 1e-3 * peak))[0]
    t_stop = h[low[0], T] if low.size else None
    r = measure(f"{path}/cells.dat", nx, ny, L_i)
    r1 = measure(f"{path}/cells_1.dat", nx, ny, L_i)
    r05 = measure(f"{path}/cells_0.5.dat", nx, ny, L_i)
    return dict(a=a, E=E_MPa, L_i=L_i, H_i=H_i, nx=nx, ny=ny, steps=int(h[-1, STEP]),
                t=h[-1, T], seconds=sec, dm=abs(h[-1, MASS] / h[0, MASS] - 1.0),
                s_empty=h[:, SEMPTY].max(), ke_end=ke[-1] / peak, t_stop=t_stop,
                runout=r["runout"], H_f=r["H_f"], slope=r["slope"],
                runout_05=r05["runout"], runout_1=r1["runout"], hc=r["hc"])


def check_parallel(work):
    """Two processors against serial: a = 0.5, 1 MPa, to t = 0.1 s, split
    across the column at y = 0.06 m (the shear band and the top layers cross
    it). Scaled differences <= 1e-10 (velocity by 1 m/s, stress by rho g
    H_i). Later the runs part: switches of the cap (bisection) and rule 9
    (a gradient comparison) flip on round-off, which the plastic flow
    amplifies (decision 0032)."""
    ps = setup(work, 0.5, 1.0, 0.1, "parSer")
    pp = setup(work, 0.5, 1.0, 0.1, "parPar", nProcs=2)
    dp = f"{pp}/system/decomposeParDict"
    txt = open(dp).read()
    open(dp, "w").write(txt.replace("n (2 1 1)", "n (1 2 1)"))
    vc.run(ps)
    vc.run(pp, 2)
    a, b = vc.cells(ps), vc.cells(pp, 2)
    sg = RHO * 9.81 * 0.1
    scale = {"alpha": (2, 1.0), "rho": (3, RHO), "Ux": (4, 1.0), "Uy": (5, 1.0),
             "s_xx": (6, sg), "s_yy": (7, sg), "s_xy": (8, sg)}
    if a.shape != b.shape:
        print("\n  2 processors: cell counts differ: FAIL")
        return False
    d = {k: np.abs(a[:, i] - b[:, i]).max() / sc for k, (i, sc) in scale.items()}
    moved = np.abs(a[:, 4]).max()
    capped = np.loadtxt(f"{pp}/history.dat")[:, NCAP].sum()
    good = max(d.values()) <= 1e-10 and moved > 0.01 and capped > 0
    print(f"\n  2 processors against serial (a = 0.5, 1 MPa, t = 0.1 s, split at y = 0.06 m, "
          f"largest |Ux| {moved:.3f} m/s, {capped:.0f} face-steps capped): largest scaled "
          "difference "
          + ", ".join(f"{k} {v:.0e}" for k, v in d.items())
          + f" (<= 1e-10): {'PASS' if good else 'FAIL'}")
    return good


def main():
    work = sys.argv[1]
    cases = [(0.5, 0.1), (0.5, 1.0), (0.5, 10.0), (2.0, 1.0)]
    if len(sys.argv) > 2:
        cases = [tuple(float(x) for x in s.split(":")) for s in sys.argv[2:]]
    t_all = time.perf_counter()
    ok = True
    print("Granular column collapse, T0 Mohr-Coulomb (c = 0, phi 35, psi 0), 10 mm cells,")
    print(f"smooth left wall, floor friction {DELTA_B:g} deg, from K0 rest, shear cap and")
    print(f"rule 9 on, neither rule 12 nor 13, t_max {T_MAX} s")
    for a, E in cases:
        r = collapse(work, a, E)
        lo, hi = BAND[a]
        in_band = lo <= r["runout"] <= hi
        top = abs(r["H_f"] - r["H_i"]) <= CELL if a < 1.0 else True
        basic = r["dm"] <= 1e-10 and r["s_empty"] == 0.0 and np.isfinite(r["runout"])
        checked = E == 1.0 or CHECK_ALL
        good = basic and (not checked or (in_band and top and r["slope"] <= 37.0))
        ok &= good
        stop = (f"stops at t = {r['t_stop']:.2f} s" if r["t_stop"] is not None
                else f"still moving at t = {r['t']:.2f} s (KE {r['ke_end']:.0e} of peak)")
        print(f"\n  a = {a:g}, E = {E:g} MPa: {r['nx']} x {r['ny']} cells, {r['steps']} steps "
              f"({r['seconds']:.0f} s); {stop}")
        print(f"    run-out {r['runout']:.2f} (at 0.5 s {r['runout_05']:.2f}, at 1 s "
              f"{r['runout_1']:.2f}; band {lo} to {hi}, prototype {PROTO.get((a, E), 'n/a')}); "
              f"H_f {r['H_f']:.4f} m (H_i {r['H_i']:.3f}); fitted slope {r['slope']:.1f} deg "
              f"(<= 37); mass change {r['dm']:.0e}; largest stress in empty cells "
              f"{r['s_empty']:.0e} Pa")
        print("    column heights (cells): "
              + " ".join(f"{x / CELL:.1f}" for x in r["hc"][:int(r['L_i'] * (1 + hi) / CELL) + 2]))
        verdict = "PASS" if good else "FAIL"
        print(f"    {'checked' if checked else 'recorded'}: {verdict}")
    if len(sys.argv) <= 2:
        ok &= check_parallel(work)
    print(f"\nColumn collapse (C++) checks: {'ALL PASS' if ok else 'FAILED'} "
          f"({time.perf_counter() - t_all:.0f} s)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
