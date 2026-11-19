"""WP 6.1 done-when check: MPI runs with Lagrangian bodies (decision 0080).
Run by tests/test_mpi_bodies.sh after celVoidCheck is built:

    python3 mpi_checks.py WORKDIR [case ...]

Each case runs serial, on N processors (the split goes through the body or
the outlet), and serial again with maxCo changed by 1 part in 1e12 (a
round-off perturbation). The explicit cycle has switches (shear and normal
caps, void thresholds, cut-cell mixing) that a round-off difference can flip,
after which runs part like two runs with different round-off (decisions
0032, 0046). So "1 and N cores agree" is checked in two ways:

- exact phase: the first STEPS_EXACT steps of history.dat (mass, momentum,
  body pose and force, ...) match serial to 1e-12 of each column's largest
  value, or to 10 times the perturbed serial run's difference where that is
  larger (a stress-free fill sits at the yield surface, so its cap counts
  and momentum follow round-off from step 2); a missing reduction or a
  wrong processor-boundary value shows here at once;
- outcome: the case result on N processors is within the tolerance of
  serial, and the perturbed serial run is printed beside it to show the
  round-off band.

Cases (10 mm cells, T0 rule set, all short):
  drop     free steel block dropped onto a T0 layer (WP 4.3d), 2 processors
           split at x = 0.15 m through the block, 0.3 s: mean force of the
           material on the block over the last 0.05 s within 2 % of serial,
           block height within 0.2 mm, mass to 1e-10, no hand-over loss.
  deform   steel footing as hex8 elements pushed into T0 (WP 5.3a), 2
           processors through the footing, push phase 0.2 s: mean force over
           0.06 to 0.19 s within 5 % of serial (the WP 5.3a
           tolerance; the perturbed serial run parts by about 2.5 %), mass to 1e-10.
  silo     T0 tank with STL walls and floor and an outlet (WP 4.5), E 0.1
           MPa, mu(I), 4 processors (2 x 2: the split runs through the
           outlet and the floor), opens at 0.3 s, runs 0.8 s more: rate within
           5 % of serial, nothing through the closed outlet (4e-14), tank
           plus discharged constant to 1e-12.

Prints PASS/FAIL per case, exits non-zero on any failure.
"""

import math
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
for d in ("silo", "void", "lagrange", "stlwalls", "freebody", "deformable", "outlet"):
    sys.path.insert(0, os.path.join(HERE, "..", d))
import cases as cs  # noqa: E402
import void_checks as vc  # noqa: E402
import stl_push_checks as sp  # noqa: E402
import free_body_checks as fb  # noqa: E402
import deformable_push_checks as dp  # noqa: E402
import outlet_checks as oc  # noqa: E402
import stl_silo_checks as ss  # noqa: E402

STEPS_EXACT = 100
EXACT_TOL = 1e-12
PERTURB = 1 + 1e-12


def decompose(path, n):
    split = {2: "(2 1 1)", 4: "(2 2 1)"}[n]
    cs.write(f"{path}/system/decomposeParDict",
             cs.HEADER.format(cls="dictionary", obj="decomposeParDict")
             + f"numberOfSubdomains {n};\nmethod simple;\ncoeffs {{ n {split}; }}\n")


def set_maxco(path, factor):
    fn = f"{path}/system/voidCheckDict"
    s = open(fn).read()
    m = re.search(r"^maxCo (.*);$", s, flags=re.M)
    s = s.replace(m.group(0), f"maxCo {float(m.group(1)) * factor!r};")
    open(fn, "w").write(s)


def exact_phase(ser, par, pert):
    """First STEPS_EXACT steps, per column scaled by the column's largest
    value over the serial run: largest difference of the parallel run, of
    the perturbed serial run, and whether every column of the parallel run
    is within max(EXACT_TOL, 10 x the perturbed run's difference)."""
    k = min(STEPS_EXACT, len(ser), len(par), len(pert))
    scale = np.abs(ser).max(0)
    scale[scale == 0] = 1.0
    dp = (np.abs(ser[:k] - par[:k]) / scale).max(0)
    dr = (np.abs(ser[:k] - pert[:k]) / scale).max(0)
    return dp.max(), dr.max(), bool(np.all(dp <= np.maximum(EXACT_TOL, 10 * dr)))


def exact_line(e):
    return (f"first {STEPS_EXACT} steps: largest scaled difference {e[0]:.1e}, perturbed "
            f"serial {e[1]:.1e} (each column <= max({EXACT_TOL}, 10 x perturbed): "
            f"{'yes' if e[2] else 'no'})")


def run_three(work, name, build, nProcs):
    """Serial, parallel and perturbed serial runs of one case; returns
    their history arrays and the wall time of each."""
    out, times = {}, {}
    for tag, n, f in (("ser", 1, 1.0), ("par", nProcs, 1.0), ("pert", 1, PERTURB)):
        p = build(f"{work}/{name}_{tag}")
        if n > 1:
            decompose(p, n)
        if f != 1.0:
            set_maxco(p, f)
        times[tag] = vc.run(p, n)
        out[tag] = (p, np.loadtxt(f"{p}/history.dat"))
    return out, times


def report(name, ok, lines):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}")
    for line in lines:
        print(f"    {line}")
    return ok


# --------------------------------------------------------------------------

def drop_case(path):
    p = vc.make_case(path, 30, 30, vc.MC + sp.T0_RULES,
                     {"fill": [-1.0, -1.0, 1.0, fb.DEPTH], "U0": [0, 0, 0],
                      "initK0": vc.NU / (1 - vc.NU), "tEnd": 0.3, "maxCo": 0.4})
    m, _ = fb.block_mass()
    fb.add_free_body(p, "block", (fb.X0, fb.Y0, -0.01), (fb.X1, fb.Y0 + fb.HEIGHT, 0.02),
                     m, 0.0, (0.0, 0.0), 0.0, 20.0)
    return p


def drop_result(h):
    m, _ = fb.block_mass()
    late = h[:, sp.T] >= h[-1, sp.T] - 0.05
    return dict(force=-h[late, sp.FYB].mean() / (m * vc.G_ACC), y=h[-1, fb.CY],
                dm=abs(h[-1, sp.MASS] / h[0, sp.MASS] - 1), lost=h[:, sp.LOST].max())


def check_drop(work):
    print("[1] Free block dropped onto T0, serial against 2 processors (split through it)")
    r, t = run_three(work, "drop", drop_case, 2)
    ex = exact_phase(r["ser"][1], r["par"][1], r["pert"][1])
    s, q, e = (drop_result(r[k][1]) for k in ("ser", "par", "pert"))
    df = abs(q["force"] / s["force"] - 1)
    dy = abs(q["y"] - s["y"])
    ok = (ex[2] and df <= 0.02 and dy <= 2e-4 and q["dm"] <= 1e-10
          and q["lost"] == 0)
    return report("drop", ok, [
        exact_line(ex),
        f"force / weight over the last 0.05 s: serial {s['force']:.4f}, 2 processors "
        f"{q['force']:.4f} ({100 * df:.2f} %, <= 2 %), perturbed serial {e['force']:.4f}",
        f"block centre height: 2 processors - serial {1000 * dy:.3f} mm (<= 0.2), "
        f"perturbed serial - serial {1000 * abs(e['y'] - s['y']):.3f} mm",
        f"2 processors: mass change {q['dm']:.0e} (<= 1e-10), lost in hand-overs "
        f"{q['lost']:.0e} (0)",
        f"wall time serial {t['ser']:.0f} s, 2 processors {t['par']:.0f} s"])


# --------------------------------------------------------------------------

def deform_case(path):
    p = vc.make_case(path, 30, 30, vc.MC + sp.T0_RULES,
                     {"fill": [-1.0, -1.0, 1.0, sp.DEPTH], "U0": [0, 0, 0],
                      "initK0": vc.NU / (1 - vc.NU), "tEnd": 0.19, "maxCo": 0.4})
    dp.add_deformable(p, 210e9, 0.3, 7850.0, dp.table(2 * sp.T_PUSH))
    return p


def deform_result(h):
    sel = (h[:, sp.T] >= dp.WINDOW[0]) & (h[:, sp.T] <= dp.WINDOW[1])
    return dict(force=h[sel, sp.FYB].mean(),
                dm=abs(h[-1, sp.MASS] / h[0, sp.MASS] - 1), lost=h[:, sp.LOST].max())


def check_deform(work):
    print("\n[2] Steel hex8 footing pushed into T0, serial against 2 processors")
    r, t = run_three(work, "deform", deform_case, 2)
    ex = exact_phase(r["ser"][1], r["par"][1], r["pert"][1])
    s, q, e = (deform_result(r[k][1]) for k in ("ser", "par", "pert"))
    df = abs(q["force"] / s["force"] - 1)
    ok = ex[2] and df <= 0.05 and q["dm"] <= 1e-10 and q["lost"] == 0
    return report("deform", ok, [
        exact_line(ex),
        f"mean force {dp.WINDOW[0]} to {dp.WINDOW[1]} s: serial {-s['force']:.2f} N, "
        f"2 processors {-q['force']:.2f} N ({100 * df:.2f} %, <= 5 %), perturbed serial "
        f"{-e['force']:.2f} N",
        f"2 processors: mass change {q['dm']:.0e} (<= 1e-10), lost {q['lost']:.0e} (0)",
        f"wall time serial {t['ser']:.0f} s, 2 processors {t['par']:.0f} s"])


# --------------------------------------------------------------------------

def silo_case(path):
    work, name = os.path.split(path)
    return ss.make(work, name, "stl", 0.4, t_dis=0.8)


def check_silo(work):
    print("\n[3] T0 tank with STL walls and an outlet, serial against 4 processors (2 x 2)")
    r, t = run_three(work, "silo", silo_case, 4)
    ex = exact_phase(r["ser"][1], r["par"][1], r["pert"][1])
    s, q, e = (ss.result(r[k][0]) for k in ("ser", "par", "pert"))
    dr = abs(q["rate"] / s["rate"] - 1)
    ok = (ex[2] and math.isfinite(dr) and dr <= 0.05 and q["leak"] <= 4e-14
          and q["cons"] <= 1e-12)
    return report("silo", ok, [
        exact_line(ex),
        f"rate: serial {s['rate']:.2f}, 4 processors {q['rate']:.2f} kg/s/m "
        f"({100 * dr:.2f} %, <= 5 %), perturbed serial {e['rate']:.2f}",
        f"4 processors: closed outlet leak {q['leak']:.0e} (<= 4e-14), tank plus "
        f"discharged {q['cons']:.0e} (<= 1e-12)",
        f"wall time serial {t['ser']:.0f} s, 4 processors {t['par']:.0f} s"])


def main():
    work = sys.argv[1]
    checks = {"drop": check_drop, "deform": check_deform, "silo": check_silo}
    only = sys.argv[2:] or list(checks)
    t0 = time.perf_counter()
    ok = True
    for name in only:
        ok = checks[name](work) and ok
    print(f"\nWP 6.1 checks: {'ALL PASS' if ok else 'FAILED'} "
          f"({time.perf_counter() - t0:.0f} s)")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
