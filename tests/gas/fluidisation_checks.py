"""WP F.7 done-when part 2: gas through and around the grains (decisions 0125, 0132).
Cases in tests/gas/fluid_cases.py and tests/gas/venting_cases.py; every run
uses the ergunWenYu drag law (Sia, 2026-11-24) and `model gasFlow` with the
pore gas momentum on (poreInertia yes, the default).

R  Pore gas momentum: a frozen bed (d 5 mm, e 0.40), the gas at rest, a
   pressure step across it: the flux follows e rho_g dq/dt = e^2 G - K q
   (tools: RK4) within 2 % from t = 0.4 tau, and its end within 0.2 %
   (backward Euler lags by about dt/(2 tau)).
O  Onset of fluidisation: sand (d 0.2 mm, e_mf 0.38, 5 cm), air ramped
   through the floor: U_mf where the pressure drop reaches the bed weight,
   within 10 % of Wen and Yu's correlation and within 2 % of the drag law's
   own U_mf; the bed rises above it (lift speed extrapolated to 0 within
   10 % of Wen-Yu); the pressure drop after lift within 10 % of the weight.
S  Dilute settling: a slab at e 0.9 (d 50 um) falls through the closed
   column's air: its speed within 10 % of Richardson and Zaki (U_t e^n,
   Dallavalle U_t, Rowe's n) and within 5 % of the drag law's own value.
C  Stage A reproduced (decision 0125 item 6): F.4's venting column and
   outlet box with the pore gas momentum against model poreGas (Darcy, void
   at p0): degree of venting and discharged mass within 1 %.
P  2 and 4 processors against serial on the onset run: 2 x 1 (cut through
   the bed and the free surface, vertical) and 1 x 4 (cuts in the bed, at
   the surface and in the void). Up to the onset: p, Ug and U within 1e-10
   of their largest |value| (pressure solve to 1e-12 Pa); gas mass less the
   inflow conserved to 1e-12 on every count. After the onset the lifting
   bed parts serial from parallel (round-off amplified, up to 1e-2:
   reported, F.9; decision 0132).
K  Cost per step with gas on against gas off (report; F.9 target 1.5).

    python3 fluidisation_checks.py WORKDIR
"""

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import gas_reference as gr  # noqa: E402
import fluid_cases as fl  # noqa: E402
import gas_flow_checks as gf  # noqa: E402
import pore_cases as pc  # noqa: E402
import venting_cases as vc  # noqa: E402

NBED = fl.NX * (fl.NB + fl.NV)
res = []


def report(name, good, detail):
    res.append(bool(good))
    print(f"  {name}: {detail}: {'PASS' if good else 'FAIL'}", flush=True)


def main(work):
    t_all = time.perf_counter()
    secs, nst = {}, {}

    # R
    print("R: pore gas momentum, flux after a pressure step (frozen bed)")
    w = f"{work}/R"
    fl.make_relax(w, end=0.02, write_interval=0.002)
    pc.run(w)
    ts = pc.times(w)
    q_ref, q_inf, tau = fl.relax_reference([float(t) for t in ts])
    err, last = 0.0, 0.0
    for t, r in zip(ts, q_ref):
        q = pc.field(w, t, "qGas", fl.N_RLX)[:, 1]
        e = float(np.abs(q / r - 1).max())
        if float(t) >= 0.4 * tau:
            err = max(err, e)
        last = e
    report("R relaxation", err <= 0.02 and last <= 2e-3,
           f"tau {tau:.4g} s, q_inf {q_inf:.4g} m/s: largest {100 * err:.2f} % from 0.4 tau (<= 2 %), "
           f"at {ts[-1]} s {100 * last:.3f} % (<= 0.2 %)")

    # O
    print("\nO: onset of fluidisation")
    w = f"{work}/O1"
    fl.make_bed(w, extra_gas="pTolerance 1e-12;\n")
    secs["O1"] = pc.run(w)
    nst["O1"] = gf.steps(w)
    u_dp, u_lift, rows = fl.onset(w)
    e_wy, e_m = u_dp / fl.U_WY - 1, u_dp / fl.U_MODEL - 1
    report("O U_mf from the pressure drop", abs(e_wy) <= 0.10 and abs(e_m) <= 0.02,
           f"{u_dp:.5f} m/s: {100 * e_wy:+.2f} % of Wen-Yu {fl.U_WY:.5f} (<= 10 %), "
           f"{100 * e_m:+.2f} % of the drag law's {fl.U_MODEL:.5f} (<= 2 %)")
    report("O bed lifts", abs(u_lift / fl.U_WY - 1) <= 0.10,
           f"rise speed extrapolated to 0 at {u_lift:.5f} m/s ({100 * (u_lift / fl.U_WY - 1):+.2f} % of "
           f"Wen-Yu); at the end it rises at {rows[-1, 3]:.4f} m/s, centre of mass up "
           f"{1e3 * (rows[-1, 4] - rows[0, 4]):.3f} mm")
    after = rows[rows[:, 1] > 1.15 * u_dp]
    dpa = float(after[:, 2].mean())
    report("O pressure drop after lift", abs(dpa / fl.W_BED - 1) <= 0.10,
           f"mean {dpa:.1f} Pa against the weight {fl.W_BED:.1f} Pa ({100 * (dpa / fl.W_BED - 1):+.2f} %)")
    report("O gas mass balance", gf.drift(w) <= 1e-12, f"{gf.drift(w):.1e}")

    # S
    print("\nS: hindered settling of a dilute slab")
    w = f"{work}/S"
    fl.make_slab(w)
    pc.run(w)
    v = np.array([fl.slab_state(w, t)[0] for t in pc.times(w)])
    vt = -float(v[-3:].mean())
    report("S settling speed", abs(vt / fl.V_RZ - 1) <= 0.10 and abs(vt / fl.V_MODEL - 1) <= 0.05,
           f"{vt:.5f} m/s (steady to {abs(v[-1] - v[-3]) / vt:.1e}): {100 * (vt / fl.V_RZ - 1):+.2f} % of "
           f"Richardson-Zaki {fl.V_RZ:.5f} (<= 10 %), {100 * (vt / fl.V_MODEL - 1):+.2f} % of the drag "
           f"law's {fl.V_MODEL:.5f} (<= 5 %)")
    report("S gas mass balance", gf.drift(w) <= 1e-12, f"{gf.drift(w):.1e}")

    # C
    print("\nC: Stage A reproduced (against model poreGas)")
    TS = (0.1, 0.5)
    col = {}
    for model in ("poreGas", "gasFlow"):
        w = f"{work}/C_{model}"
        vc.make_column(w, end=vc.t_of(0.5), write_interval=vc.t_of(0.1))
        gf.set_model(w, model)
        pc.run(w)
        col[model] = gf.venting(w, TS)
    for T in TS:
        Uf, Up = col["gasFlow"][T], col["poreGas"][T]
        report(f"C venting T = {T}", abs(Uf - Up) <= 0.01 * Up,
               f"U {Uf:.5f} against {Up:.5f} ({100 * (Uf / Up - 1):+.3f} %, <= 1 %)")
    box = {}
    for model in ("poreGas", "gasFlow", "off"):
        w = f"{work}/C_box_{model}"
        vc.make_outlet(w, end=0.15, write_interval=0.05, gas=model != "off")
        if model != "off":
            gf.set_model(w, model)
        secs["B" + model] = pc.run(w, set_fields=True)
        nst["B" + model] = gf.steps(w)
        box[model] = w
    mf, mp = gf.discharged(box["gasFlow"]), gf.discharged(box["poreGas"])
    report("C outlet box discharge", mf > 0 and abs(mf - mp) <= 0.01 * mp,
           f"{mf:.6g} kg against {mp:.6g} kg ({100 * (mf / mp - 1):+.3f} %, <= 1 %); gas mass "
           f"{gf.drift(box['gasFlow']):.1e}")

    # P
    print("\nP: onset run on 2 and 4 processors against serial")
    out = {1: f"{work}/O1"}
    t_on = u_dp / fl.U_END * fl.T_RAMP
    for n, dec in ((2, (2, 1, 1)), (4, (1, 4, 1))):
        w = f"{work}/O{n}"
        fl.make_bed(w, decomp=dec, extra_gas="pTolerance 1e-12;\n")
        secs[f"O{n}"] = pc.run(w, n)
        out[n] = w
        b, per, aft = fl.compare_runs(out[1], w, NBED, t_on)
        report(f"P {n} processors {dec[0]} x {dec[1]}, up to the onset ({t_on:.3f} s)",
               max(b.values()) <= 1e-10 and gf.drift(w) <= 1e-12,
               f"p within {b['p']:.1e}, Ug {b['Ug']:.1e}, U {b['U']:.1e} of their largest (each "
               f"write's own largest: {per['p']:.1e}, {per['Ug']:.1e}, {per['U']:.1e}); gas mass "
               f"{gf.drift(w):.1e}")
        print(f"    after the onset (report, F.9): p {aft['p']:.1e}, Ug {aft['Ug']:.1e}, U {aft['U']:.1e}")
    print(f"  wall times: {secs['O1']:.1f} / {secs['O2']:.1f} / {secs['O4']:.1f} s on 1 / 2 / 4 "
          f"processors ({NBED} cells)")

    # K
    print("\nK: cost per step, gas on against gas off (report; F.9 target 1.5)")
    w = f"{work}/K_off"
    fl.make_bed(w, gas=False)
    secs["Koff"] = pc.run(w)
    nst["Koff"] = gf.steps(w)

    def per(k):
        return secs[k] / nst[k]

    print(f"  fluidised bed ({NBED} cells; dry, the bed only sits): {per('O1') / per('Koff'):.2f}x")
    print(f"  outlet box ({vc.OUT_NX * vc.OUT_NY} cells): poreGas {per('BporeGas') / per('Boff'):.2f}x, "
          f"gasFlow with pore momentum {per('BgasFlow') / per('Boff'):.2f}x")

    print(f"\ntest_fluidisation: {'PASS' if all(res) else 'FAIL'} ({sum(res)}/{len(res)}), "
          f"{time.perf_counter() - t_all:.0f} s")
    return 0 if all(res) else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
