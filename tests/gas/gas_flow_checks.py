"""WP F.6 done-when: gas momentum in the void (decisions 0125, 0131).
Cases in tests/gas/flow_cases.py and tests/gas/venting_cases.py; every run
uses the ergunWenYu drag law and `model gasFlow` (void cells carry a gas
velocity Ug and their own pressure).

A  Plane channel, void only, from rest to steady: Ug against Poiseuille,
   largest error within 1 % of u_max and mean within 1 %; the time step
   equals the solid's step without gas, step by step; gas mass less the
   inflow constant to 1e-12.
B  Still air in a closed box of non-orthogonal cells: |Ug| below 1e-10 m/s.
C  Bed and void together (Stage A reproduced, decision 0125 item 6):
   venting column of WP F.4 with the void now solved: within 1 % of
   Terzaghi and of the poreGas run (void at p0); Mohr-Coulomb outlet box
   (material and gas leaving): discharged mass within 1 % of poreGas;
   gas mass balance 1e-12; time step unchanged.
P  2 and 4 processors against serial: channel (cuts along and across it,
   through the inlet and outlet), venting column (cuts in the bed and at
   the free surface), Mohr-Coulomb outlet box 2 x 1 and 2 x 2 (cuts
   through the outlet, the bed and the free surface, material leaving).
   p within 1e-10 of the largest |p| (pressure solve to 1e-12 of it),
   Ug and U within 1e-10 of their largest; gas mass 1e-12 on every count.
K  Cost per step, gas on (gasFlow and poreGas) against gas off (report;
   F.9 sets the 1.5 limit).

    python3 gas_flow_checks.py WORKDIR
"""

import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import gas_reference as gr  # noqa: E402
import flow_cases as fc  # noqa: E402
import pore_cases as pc  # noqa: E402
import venting_cases as vc  # noqa: E402

NCH = fc.NX * fc.NY
NCOL = vc.NB + vc.NV
NOUT = vc.OUT_NX * vc.OUT_NY
res = []


def report(name, good, detail):
    res.append(bool(good))
    print(f"  {name}: {detail}: {'PASS' if good else 'FAIL'}", flush=True)


def set_model(path, model):
    f = f"{path}/constant/gasProperties"
    txt = open(f).read().replace("model   poreGas;", f"model   {model};")
    open(f, "w").write(txt)


def drift(path):
    m, ni, _, _ = pc.mass_log(path)
    return float(np.abs((m - ni) - (m[0] - ni[0])).max() / m[0])


def log(path):
    return open(f"{path}/log.celFoam").read()


def dts(path):
    return re.findall(r"^Time = \S+  deltaT = (\S+)", log(path), re.M)


def steps(path):
    return len(dts(path))


def discharged(path):
    return float(re.findall(r"Discharged mass = (\S+)", log(path))[-1])


def vfield(path, t, name, n):
    """pc.field, also for uniform vector fields"""
    txt = open(f"{path}/{t}/{name}").read()
    m = re.search(r"internalField\s+uniform\s+\((\S+) (\S+) (\S+)\);", txt)
    if m:
        return np.tile([float(x) for x in m.groups()], (n, 1))
    return pc.field(path, t, name, n)


def rel(a, b):
    return float(np.abs(a - b).max() / max(np.abs(a).max(), 1e-300))


def compare(a, b, n, fields):
    out = {f: 0.0 for f in fields}
    for t in pc.times(a):
        for f in fields:
            out[f] = max(out[f], rel(vfield(a, t, f, n), vfield(b, t, f, n)))
    return out


def venting(path, Ts):
    c = vc.c_of(vc.E_MOD)
    out = {}
    for t in pc.times(path):
        T = float(t) * c / vc.HB ** 2
        for x in Ts:
            if abs(T - x) < 1e-6:
                out[x] = vc.degree(path, t)
    return out


def main(work):
    print(f"Channel: H {fc.H} m, L {fc.L} m, {fc.NY} x {fc.NX} cells, DP {fc.DP} Pa, "
          f"u_max {fc.U_MAX:.5g} m/s, Re {fc.U_MEAN * fc.H / fc.NU_G:.3g}, to {fc.END:.4g} s "
          f"(1.2 H^2/nu)")
    secs, nst = {}, {}

    # A
    print("A: plane channel against Poiseuille")
    a = f"{work}/A"
    # Two pressure corrections per step (PISO; the default is one)
    fc.make_channel(a, extra_gas="pTolerance 1e-15;\nnCorr 2;\n")
    secs["A1"] = pc.run(a)
    nst["A1"] = steps(a)
    C = pc.cell_centres(a)
    ue = fc.u_exact(C[:, 1])
    t = pc.times(a)[-1]
    U = pc.field(a, t, "Ug", NCH)
    emax = float(np.abs(U[:, 0] - ue).max()) / fc.U_MAX
    emean = float(U[:, 0].mean()) / fc.U_MEAN - 1
    report("Poiseuille profile", emax <= 0.01, f"largest error {emax:.3%} of u_max")
    report("Poiseuille mean velocity", abs(emean) <= 0.01, f"{emean:+.3%}")
    report("no cross flow", float(np.abs(U[:, 1]).max()) <= 1e-6 * fc.U_MAX,
           f"largest |v| {np.abs(U[:, 1]).max():.1e} m/s")
    report("gas mass balance (A)", drift(a) <= 1e-12, f"{drift(a):.1e}")
    a0 = f"{work}/A_off"
    fc.make_channel(a0, gas=False)
    secs["A0"] = pc.run(a0)
    nst["A0"] = steps(a0)
    report("time step of the solid (channel)", dts(a) == dts(a0),
           f"{len(dts(a))} steps, deltaT {dts(a)[0]} s, same as without gas")

    # B
    print("B: still air, closed box of non-orthogonal cells")
    b = f"{work}/B"
    fc.make_box(b, end=0.05, write_interval=0.025)
    pc.run(b)
    um = max(float(np.abs(vfield(b, t, "Ug", NCH)).max()) for t in pc.times(b))
    report("still air stays still", um < 1e-10, f"largest |Ug| {um:.1e} m/s")

    # C
    print("C: bed and void together (Stage A reproduced)")
    TS = (0.05, 0.1, 0.2, 0.5)
    col = {}
    for model in ("poreGas", "gasFlow"):
        w = f"{work}/C_{model}"
        vc.make_column(w, end=vc.t_of(0.5), write_interval=vc.t_of(0.05))
        set_model(w, model)
        secs["Cc" + model] = pc.run(w)
        nst["Cc" + model] = steps(w)
        col[model] = venting(w, TS)
        cf = w
    w = f"{work}/C_off"
    vc.make_column(w, end=vc.t_of(0.5), write_interval=vc.t_of(0.05), gas=False)
    secs["Cc0"] = pc.run(w)
    nst["Cc0"] = steps(w)
    for T in TS:
        Ut, Uf, Up = gr.terzaghi_degree(T), col["gasFlow"][T], col["poreGas"][T]
        report(f"venting T = {T}", abs(Uf - Ut) <= 0.01 * Ut and abs(Uf - Up) <= 0.01 * Up,
               f"U {Uf:.4f}: {(Uf - Ut) / Ut:+.2%} of Terzaghi, {(Uf - Up) / Up:+.2%} of poreGas")
    report("gas mass balance (column)", drift(cf) <= 1e-12, f"{drift(cf):.1e}")
    # Report only: the column's step is the solid's, but its surface cell's
    # density moves with the gas pressure, so the stable step differs from
    # the dry run by about 1e-6, and by up to one step in a few hundred where
    # a write time shortens it (the exact check is A and the outlet box)
    d1, d0 = np.array(dts(cf), float), np.array(dts(w), float)
    if len(d1) == len(d0):
        dd = np.abs(d1 - d0) / d0
        print(f"  column time step (report): {len(d1)} steps as without gas, median "
              f"difference {np.median(dd):.1e}, largest {dd.max():.1e}")
    else:
        print(f"  column time step (report): {len(d1)} steps against {len(d0)} without gas")
    um = max(float(np.abs(vfield(cf, t, "Ug", NCOL)[vc.NB:]).max()) for t in pc.times(cf))
    print(f"  largest |Ug| in the void of the column {um:.3g} m/s")
    box = {}
    for model in ("poreGas", "gasFlow", "off"):
        w = f"{work}/C_box_{model}"
        vc.make_outlet(w, end=0.15, write_interval=0.05, gas=model != "off")
        if model != "off":
            set_model(w, model)
        secs["Cb" + model] = pc.run(w, set_fields=True)
        nst["Cb" + model] = steps(w)
        box[model] = w
    mf, mp = discharged(box["gasFlow"]), discharged(box["poreGas"])
    report("outlet box discharge", mf > 0 and abs(mf - mp) <= 0.01 * mp,
           f"{mf:.6g} kg against {mp:.6g} kg with poreGas ({(mf - mp) / mp:+.2%})")
    report("gas mass balance (outlet box)", drift(box["gasFlow"]) <= 1e-12,
           f"{drift(box['gasFlow']):.1e}")
    report("time step of the solid (outlet box)", dts(box["gasFlow"]) == dts(box["off"]),
           f"{len(dts(box['gasFlow']))} steps")

    # P
    print("P: 2 and 4 processors against serial")
    tol = "pTolerance 1e-12;\n"
    ch, cl, mc = {}, {}, {}
    for n, dch, dcl, dmc in ((1, (1, 1, 1), (1, 1, 1), (1, 1, 1)),
                             (2, (2, 1, 1), (1, 2, 1), (2, 1, 1)),
                             (4, (2, 2, 1), (1, 4, 1), (2, 2, 1))):
        w = f"{work}/Pch{n}"
        fc.make_channel(w, end=0.08, write_interval=0.04, decomp=dch)
        secs[f"Pch{n}"] = pc.run(w, n)
        ch[n] = w
        w = f"{work}/Pcol{n}"
        vc.make_column(w, end=vc.t_of(0.1), write_interval=vc.t_of(0.05), decomp=dcl)
        set_model(w, "gasFlow")
        # The bed is nearly at rest (largest U 4e-4 m/s): its velocity
        # agrees to 1e-10 only with the pressure to 1e-16 of p_i
        with open(f"{w}/constant/gasProperties", "a") as f:
            f.write("pTolerance 1e-14;\n")
        pc.run(w, n)
        cl[n] = w
        w = f"{work}/Pmc{n}"
        vc.make_outlet(w, end=0.1, write_interval=0.05, decomp=dmc, extra_gas=tol)
        set_model(w, "gasFlow")
        secs[f"Pmc{n}"] = pc.run(w, n, set_fields=True)
        mc[n] = w
        d = max(drift(ch[n]), drift(cl[n]), drift(mc[n]))
        report(f"gas mass balance, {n} processor(s)", d <= 1e-12,
               f"channel {drift(ch[n]):.1e}, column {drift(cl[n]):.1e}, outlet box {drift(mc[n]):.1e}")
    for n in (2, 4):
        r = compare(ch[1], ch[n], NCH, ("p", "Ug"))
        report(f"channel, {n} processors (cuts along and across)",
               max(r.values()) <= 1e-10, f"p within {r['p']:.1e}, Ug within {r['Ug']:.1e}")
        r = compare(cl[1], cl[n], NCOL, ("p", "Ug", "U"))
        report(f"column, {n} processors (cuts in the bed and at the surface)",
               max(r.values()) <= 1e-10,
               f"p within {r['p']:.1e}, Ug within {r['Ug']:.1e}, U within {r['U']:.1e}")
        r = compare(mc[1], mc[n], NOUT, ("p", "Ug", "U"))
        report(f"Mohr-Coulomb outlet box, {n} processors (outlet, bed, surface)",
               max(r.values()) <= 1e-10,
               f"p within {r['p']:.1e}, Ug within {r['Ug']:.1e}, U within {r['U']:.1e}")
    print(f"  wall times, channel 0.08 s: {secs['Pch1']:.2f} / {secs['Pch2']:.2f} / "
          f"{secs['Pch4']:.2f} s; outlet box 0.1 s: {secs['Pmc1']:.2f} / {secs['Pmc2']:.2f} / "
          f"{secs['Pmc4']:.2f} s on 1 / 2 / 4 processors ({NCH} and {NOUT} cells: too small to scale)")

    # K
    print("K: cost per step, gas on against gas off (report; F.9 target 1.5)")

    def per(k):
        return secs[k] / nst[k]

    print(f"  venting column ({NCOL} cells): poreGas {per('CcporeGas') / per('Cc0'):.2f}x, "
          f"gasFlow {per('CcgasFlow') / per('Cc0'):.2f}x")
    print(f"  outlet box ({NOUT} cells): poreGas {per('CbporeGas') / per('Cboff'):.2f}x, "
          f"gasFlow {per('CbgasFlow') / per('Cboff'):.2f}x")
    print(f"  channel ({NCH} cells, no material: the solid step costs almost nothing): "
          f"{per('A1') / per('A0'):.1f}x")

    n_fail = res.count(False)
    print(f"{len(res) - n_fail} of {len(res)} checks passed")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
