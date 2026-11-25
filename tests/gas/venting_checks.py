"""WP F.4 done-when: gas with a free surface, void and an outlet (decisions
0125, 0129). Cases in tests/gas/venting_cases.py; every run uses the
ergunWenYu drag law, two-way coupling and a freeSurface solid.

A  Venting column: fine powder bed 0.5 m (21 cells) under 7 void cells,
   gas at p_i = 100 Pa over p0 in the bed, the void and the top at p0.
   Degree of venting U = 1 - mean(p)/p_i against Terzaghi at T = 0.05,
   0.1, 0.2 and 0.5 within 3 %, c = k/(mu (e/p0 + 1/M)); the void cells
   hold p0 exactly; gas mass less the inflow (boundaries and void) constant
   to 1e-12. Also reported: the soft bed of WP F.3 (E 1 MPa), whose surface
   cell creeps (v1 partly filled cell).
B  Outlet box (Mohr-Coulomb, 1 mm grains): the outlet is sealed for the gas
   until it opens (no gas crosses the boundaries before), holds p0 at every
   write, and gas crosses it after; gas mass balance 1e-12; discharged mass
   with gas within 1 % of the run without gas (coarse grains: the gas barely
   matters).
P  2 and 4 processors against serial. Column to T = 0.1, cuts across it
   in the bed and at the free surface. Elastic outlet box (outlet open from
   0.02 s, run to 0.05 s), cuts through the outlet
   and the bed. Mohr-Coulomb outlet box to 0.05 s (material leaving) on 4
   processors (2 x 2). p within 1e-10 (of p_i, or of the largest |p|) with
   the solve converged to 1e-12 Pa, U within 1e-10, gas mass 1e-12 on
   every count.
C  Cost per step, gas on against gas off: column, outlet box and the outlet
   box at 2.5 mm cells (report; F.9 sets the 1.5 limit).

    python3 venting_checks.py WORKDIR
"""

import os
import re
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import gas_reference as gr  # noqa: E402
import pore_cases as pc  # noqa: E402
import venting_cases as vc  # noqa: E402

TS = (0.05, 0.1, 0.2, 0.5)
NCOL = vc.NB + vc.NV
NOUT = vc.OUT_NX * vc.OUT_NY
res = []


def report(name, good, detail):
    res.append(bool(good))
    print(f"  {name}: {detail}: {'PASS' if good else 'FAIL'}", flush=True)


def drift(path):
    m, ni, _, _ = pc.mass_log(path)
    return float(np.abs((m - ni) - (m[0] - ni[0])).max() / m[0])


def boundary_log(path):
    """Per step: time, gas that came in through the boundaries so far
    (net in less the void's share) [kg]."""
    txt = open(f"{path}/log.celFoam").read()
    t = [float(x) for x in re.findall(r"^Time = (\S+)", txt, re.M)]
    v = [float(x) for x in re.findall(r"^Pore gas void: \S+ cells, from the void (\S+) kg", txt, re.M)]
    _, ni, _, _ = pc.mass_log(path)
    n = min(len(t), len(v), len(ni))
    return np.array(t[:n]), ni[:n] - np.array(v[:n])


def steps(path):
    return len(re.findall(r"^Time = ", open(f"{path}/log.celFoam").read(), re.M))


def discharged(path):
    return float(re.findall(r"Discharged mass = (\S+)", open(f"{path}/log.celFoam").read())[-1])


def venting(path, E, Ts):
    c = vc.c_of(E)
    out = {}
    for t in pc.times(path):
        T = float(t) * c / vc.HB ** 2
        for x in Ts:
            if abs(T - x) < 1e-6:
                out[x] = vc.degree(path, t)
    return out


def compare(a, b, n, scale_p):
    dp = du = 0.0
    for t in pc.times(a):
        dp = max(dp, float(np.abs(pc.field(a, t, "p", n) - pc.field(b, t, "p", n)).max()) / scale_p)
        ua, ub = pc.field(a, t, "U", n), pc.field(b, t, "U", n)
        du = max(du, float(np.abs(ua - ub).max() / max(np.abs(ua).max(), 1e-300)))
    return dp, du


def main(work):
    print(f"Column: M {vc.M_MOD:.6g} Pa, c {vc.C_V:.6g} m^2/s, T = 1 at {vc.t_of(1):.4g} s, "
          f"p_i {vc.P_I} Pa (top-cell weight {vc.RHO * vc.G_ACC * vc.DY:.0f} Pa)")
    secs, nsteps = {}, {}

    # A
    print("A: venting column against Terzaghi")
    a = f"{work}/A"
    vc.make_column(a, end=vc.t_of(0.5), write_interval=vc.t_of(0.05))
    secs["A1"] = pc.run(a)
    nsteps["A1"] = steps(a)
    U = venting(a, vc.E_MOD, TS)
    for T in TS:
        Ut = gr.terzaghi_degree(T)
        report(f"T = {T}", abs(U[T] - Ut) <= 0.03 * Ut,
               f"U {U[T]:.4f} against {Ut:.4f} ({(U[T] - Ut) / Ut:+.2%})")
    pv = max(float(np.abs(pc.field(a, t, "p", NCOL)[vc.NB:]).max()) for t in pc.times(a))
    report("void cells at p0", pv == 0.0, f"largest |p| in the void {pv:.1e} Pa")
    report("gas mass balance (A)", drift(a) <= 1e-12, f"{drift(a):.1e}")
    off = f"{work}/A_off"
    vc.make_column(off, end=vc.t_of(0.5), write_interval=vc.t_of(0.5), gas=False)
    secs["A0"] = pc.run(off)
    nsteps["A0"] = steps(off)

    soft = f"{work}/A_soft"
    vc.make_column(soft, end=vc.t_of(0.5, 1e6), write_interval=vc.t_of(0.05, 1e6), E=1e6)
    pc.run(soft)
    Us = venting(soft, 1e6, TS)
    print("  soft bed (E 1 MPa, report): " + ", ".join(
        f"T {T}: {(Us[T] - gr.terzaghi_degree(T)) / gr.terzaghi_degree(T):+.2%}" for T in TS)
        + f"; surface cell alpha {pc.field(soft, pc.times(soft)[-1], 'alpha', NCOL)[vc.NB - 1]:.6f}")

    # B
    print("B: outlet box (Mohr-Coulomb, 1 mm grains, opens at "
          f"{vc.OPEN_TIME} s)")
    b = f"{work}/B"
    vc.make_outlet(b, end=0.15, write_interval=0.05)
    secs["B1"] = pc.run(b, set_fields=True)
    nsteps["B1"] = steps(b)
    t, inb = boundary_log(b)
    before = float(np.abs(inb[t < vc.OPEN_TIME - 1e-12]).max())
    after = float(np.abs(inb[-1]))
    report("outlet sealed until it opens", before == 0.0 and after > 0,
           f"gas through the boundaries before {before:.1e} kg, by the end {inb[-1]:+.3e} kg")
    pout = max(float(np.abs(pc.patch_value(b, t, "p", "outlet")).max()) for t in pc.times(b))
    report("outlet holds p0", pout == 0.0, f"largest |p| on the outlet {pout:.1e} Pa")
    report("gas mass balance (B)", drift(b) <= 1e-12, f"{drift(b):.1e}")
    boff = f"{work}/B_off"
    vc.make_outlet(boff, end=0.15, write_interval=0.15, gas=False)
    secs["B0"] = pc.run(boff, set_fields=True)
    nsteps["B0"] = steps(boff)
    mo, mf = discharged(b), discharged(boff)
    report("discharge with gas (coarse grains)", mo > 0 and abs(mo - mf) <= 0.01 * mf,
           f"{mo:.6g} kg against {mf:.6g} kg without gas ({(mo - mf) / mf:+.2%})")

    # P
    print("P: 2 and 4 processors against serial")
    col, box = {}, {}
    tol = "pTolerance 1e-12;\n"
    for n, dc, db in ((1, (1, 1, 1), (1, 1, 1)), (2, (1, 2, 1), (2, 1, 1)), (4, (1, 4, 1), (2, 2, 1))):
        w = f"{work}/Pc{n}"
        vc.make_column(w, end=vc.t_of(0.1), write_interval=vc.t_of(0.05), decomp=dc)
        secs[f"Pc{n}"] = pc.run(w, n)
        col[n] = w
        w = f"{work}/Pb{n}"
        vc.make_outlet(w, end=0.05, write_interval=0.01, decomp=db, material="linearElastic",
                       extra_gas=tol)
        pc.run(w, n, set_fields=True)
        box[n] = w
        report(f"gas mass balance, {n} processor(s)", max(drift(col[n]), drift(box[n])) <= 1e-12,
               f"column {drift(col[n]):.1e}, elastic outlet box {drift(box[n]):.1e}")
    pbox = max(float(np.abs(pc.field(box[1], t, "p", NOUT)).max()) for t in pc.times(box[1]))
    for n in (2, 4):
        dp, du = compare(col[1], col[n], NCOL, vc.P_I)
        report(f"column, {n} processors (cuts in the bed and at the surface)",
               dp <= 1e-10 and du <= 1e-10, f"p within {dp:.1e} of p_i, U within {du:.1e}")
        dp, du = compare(box[1], box[n], NOUT, pbox)
        report(f"elastic outlet box, {n} processors (cuts through the outlet)",
               dp <= 1e-10 and du <= 1e-10, f"p within {dp:.1e} of the largest |p|, U within {du:.1e}")
    print(f"  wall times, column: 1 proc {secs['Pc1']:.2f} s, 2 procs {secs['Pc2']:.2f} s, "
          f"4 procs {secs['Pc4']:.2f} s ({NCOL} cells: too small to scale)")
    # Mohr-Coulomb box with material leaving through the outlet
    mc = {}
    for gas in (False, True):
        for n, d in ((1, (1, 1, 1)), (4, (2, 2, 1))):
            w = f"{work}/Pm{int(gas)}{n}"
            vc.make_outlet(w, end=0.05, write_interval=0.05, decomp=d, gas=gas, extra_gas=tol)
            pc.run(w, n, set_fields=True)
            mc[gas, n] = w
    t = pc.times(mc[False, 1])[-1]
    u1, u4 = pc.field(mc[False, 1], t, "U", NOUT), pc.field(mc[False, 4], t, "U", NOUT)
    _, du = compare(mc[True, 1], mc[True, 4], NOUT, 1.0)
    dpm = float(np.abs(pc.field(mc[True, 1], t, "p", NOUT) - pc.field(mc[True, 4], t, "p", NOUT)).max()
                / np.abs(pc.field(mc[True, 1], t, "p", NOUT)).max())
    report(f"Mohr-Coulomb box at {t} s (material leaving), 4 processors",
           dpm <= 1e-10 and du <= 1e-10,
           f"gas on: p within {dpm:.1e}, U within {du:.1e} (gas off: U within "
           f"{np.abs(u1 - u4).max() / np.abs(u1).max():.1e})")
    report("gas mass balance, Mohr-Coulomb box on 4 processors", drift(mc[True, 4]) <= 1e-12,
           f"{drift(mc[True, 4]):.1e}")

    # C
    print("C: cost per step, gas on against gas off (report)")

    def cost(on, off):
        return (secs[on] / nsteps[on]) / (secs[off] / nsteps[off])

    print(f"  column ({NCOL} cells): {secs['A1']:.2f} s / {nsteps['A1']} steps against "
          f"{secs['A0']:.2f} s / {nsteps['A0']}: {cost('A1', 'A0'):.2f}x")
    _, _, sv, _ = pc.mass_log(b)
    print(f"  outlet box ({NOUT} cells): {secs['B1']:.2f} s / {nsteps['B1']} against "
          f"{secs['B0']:.2f} s / {nsteps['B0']}: {cost('B1', 'B0'):.2f}x, {sv.mean():.2f} solves per step")
    for gas in (False, True):
        w = f"{work}/C{int(gas)}"
        vc.make_outlet(w, end=0.04, write_interval=0.04, gas=gas, refine=4)
        secs[f"C{int(gas)}"] = pc.run(w, set_fields=True)
        nsteps[f"C{int(gas)}"] = steps(w)
    _, _, sv, _ = pc.mass_log(f"{work}/C1")
    print(f"  outlet box at 2.5 mm ({16 * NOUT} cells, to 0.04 s): {secs['C1']:.2f} s / "
          f"{nsteps['C1']} against {secs['C0']:.2f} s / {nsteps['C0']}: {cost('C1', 'C0'):.2f}x, "
          f"{sv.mean():.2f} solves per step")

    print(f"{sum(res)}/{len(res)} checks passed")
    return all(res)


if __name__ == "__main__":
    work = sys.argv[1]
    os.makedirs(work, exist_ok=True)
    for d in os.listdir(work):
        shutil.rmtree(os.path.join(work, d), ignore_errors=True)
    sys.exit(0 if main(work) else 1)
