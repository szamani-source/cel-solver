"""WP F.8 done-when: initially fluidised dam break against Roche et al.
(2008) (cases in collapse_cases.py; decision 0133). Run by
tests/manual/test_fluidised_collapse.sh, about 8 minutes on 4 cores.
Result 2026-11-25 (decision 0133, F.8b): A fails (front -67 %, run-out
-62 %; gas mass 1e-14 passes), B passes, C fails (1e-10 to 0.36 / 0.46 s,
1.6e-8 / 1.1e-8 at 1.3 s), K 2.3x.

A  serial gas run: front speed over Roche's constant-velocity phase within
   20 % of sqrt(2 g h0); run-out from the gate within 20 % of 5.75 h0;
   gas mass conserved to 1e-12 (gas mass less net inflow).
B  dry collapse, gas off (no gasProperties), plain Mohr-Coulomb (version 1's
   collapse material, no cap): the front history and the final column heights
   equal those of main before version 2 (v1.1.0, stored in
   collapse_dry_reference.txt) to 1e-12. The dry run with the cap (as the gas
   run) is reported (F.8b changed the cap's p_c rule, decision 0133).
C  2 and 4 processors (manual cuts: x = 0.1 m through the column, its free
   surface and the void; 4: also y = 0.1 m through the bed) against serial:
   p, Ug and U within 1e-10 of their largest values at every write (pressure
   solve converged to pTolerance 1e-12 Pa); gas mass to 1e-12 on each count.
   The first write where they part is reported.
K  cost per step, gas on against gas off (report; F.9 target 1.5x).

Every gas run: ergunWenYu, model gasFlow, surfaceCarry yes (decision 0133),
the cap of collapse_cases.CAP.

    python3 collapse_checks.py WORKDIR [--only A,B,C,K]
"""

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import collapse_cases as cc  # noqa: E402
import gas_flow_checks as gf  # noqa: E402
import pore_cases as pc  # noqa: E402

REF = os.path.join(HERE, "collapse_dry_reference.txt")
res = []


def report(name, good, detail):
    res.append(good)
    print(f"  {'PASS' if good else 'FAIL'} {name}: {detail}")


def write_manual(path, nProcs):
    """constant/cellDecomposition: 2 procs cut at x = 0.1 m; 4 procs also
    at y = 0.1 m."""
    rows = []
    for j in range(cc.NY):
        for i in range(cc.NX):
            x, y = (i + 0.5) * cc.CELL, (j + 0.5) * cc.CELL
            k = int(x >= 0.1)
            if nProcs == 4:
                k += 2 * int(y >= 0.1)
            rows.append(str(k))
    pc.write(f"{path}/constant/cellDecomposition",
             pc.HEADER.format(cls="labelList", obj="cellDecomposition")
             + f"{cc.N}\n(\n" + "\n".join(rows) + "\n)\n")
    pc.write(f"{path}/system/decomposeParDict",
             pc.HEADER.format(cls="dictionary", obj="decomposeParDict")
             + f"numberOfSubdomains {nProcs};\nmethod manual;\n"
             + "coeffs { dataFile \"cellDecomposition\"; }\n")


def main(work, only):
    t_all = time.perf_counter()
    os.makedirs(work, exist_ok=True)
    print(cc.describe())
    secs, nst = {}, {}
    gas_extra = "pTolerance 1e-12;\n"

    w1 = f"{work}/gas1"
    if not os.path.isdir(w1) or not os.path.exists(f"{w1}/log.wallTime"):
        cc.make(w1, extra_gas=gas_extra)
    secs["gas"] = cc.run(w1)
    nst["gas"] = cc.steps(w1)
    h = cc.history(w1)

    if "A" in only:
        print("\nA: fluidised collapse against Roche et al. (2008)")
        u = cc.front_speed(h)
        ro = (h[-1, 1] - cc.X0) / cc.H0
        drift = cc.gas_mass_drift(w1)
        print("  front x [m] at t/sqrt(h0/g) = 1, 2, ..., 9: "
              + ", ".join(f"{np.interp(k * cc.T_SCALE, h[:, 0], h[:, 1]):.3f}" for k in range(1, 10)))
        report("front speed", abs(u / cc.U_FRONT - 1) <= 0.2,
               f"{u:.3f} m/s against sqrt(2 g h0) = {cc.U_FRONT:.3f} m/s ({100 * (u / cc.U_FRONT - 1):+.1f} %)")
        report("run-out", abs(ro / cc.RUNOUT - 1) <= 0.2,
               f"{ro:.2f} h0 from the gate against 5.5 to 6 ({100 * (ro / cc.RUNOUT - 1):+.1f} % of 5.75); "
               f"stops at {cc.stop_time(h) / cc.T_SCALE:.1f} sqrt(h0/g) (Roche 6.5)")
        report("gas mass", drift <= 1e-12, f"{drift:.1e} of the gas mass")

    wd = f"{work}/dry"
    if not os.path.isdir(wd) or not os.path.exists(f"{wd}/log.wallTime"):
        cc.make(wd, gas=False)
    wv1 = f"{work}/dryV1"
    if "B" in only:
        if not os.path.isdir(wv1) or not os.path.exists(f"{wv1}/log.wallTime"):
            cc.make(wv1, gas=False, cap=False)
        cc.run(wv1)
    secs["dry"] = cc.run(wd)
    nst["dry"] = cc.steps(wd)
    if "B" in only:
        print("\nB: dry collapse, gas off, against main before version 2")
        hcap = cc.history(wd)
        print(f"  dry with the cap (report): run-out {(hcap[-1, 1] - cc.X0) / cc.H0:.2f} h0, "
              f"front speed {cc.front_speed(hcap):.3f} m/s")
        hd = cc.history(wv1)
        hc, _ = cc.columns(wv1, pc.times(wv1)[-1])
        ro = (hd[-1, 1] - cc.X0) / cc.H0
        if os.path.exists(REF):
            ref = np.loadtxt(REF)
            nh = len(hd)
            dfront = float(np.abs(ref[:nh, 1] - hd[:, 1]).max()) if len(ref) >= nh + cc.NX else 1.0
            dh = float(np.abs(ref[nh:nh + cc.NX, 1] - hc).max() / hc.max()) if len(ref) >= nh + cc.NX else 1.0
            report("dry collapse unchanged", dfront <= 1e-12 and dh <= 1e-12,
                   f"front history within {dfront:.1e} m, final heights within {dh:.1e}; "
                   f"run-out {ro:.2f} h0, front speed {cc.front_speed(hd):.3f} m/s")
        else:
            np.savetxt(REF, np.vstack([hd, np.column_stack([np.arange(cc.NX), hc])]),
                       header="WP F.8 dry dam break on main (v1.1.0): rows (t, front x) at each "
                              "write, then (column, final height)")
            report("dry collapse reference written (run on main)", False, REF)

    if "C" in only:
        print("\nC: 2 and 4 processors against serial")
        for n in (2, 4):
            w = f"{work}/gas{n}"
            if not os.path.isdir(w) or not os.path.exists(f"{w}/log.wallTime"):
                cc.make(w, extra_gas=gas_extra)
                write_manual(w, n)
            try:
                secs[f"gas{n}"] = cc.run(w, n)
            except RuntimeError as err:
                # a parallel run that parts can blow up (decision 0133):
                # compare what it wrote
                secs[f"gas{n}"] = float("nan")
                print(f"  {n} processors stopped: {str(err).splitlines()[0]}")
                if os.path.isdir(f"{w}/processor0") and pc.times(f"{w}/processor0"):
                    pc.sh(w, ["reconstructPar"], "reconstructPar")
            worst, first = 0.0, None
            done = set(pc.times(w)) if os.path.isdir(w) else set()
            if len(done) < len(pc.times(w1)):
                worst = float("inf")
            for t in [t for t in pc.times(w1) if t in done]:
                d = max(gf.rel(gf.vfield(w1, t, f, cc.N), gf.vfield(w, t, f, cc.N))
                        for f in ("p", "Ug", "U"))
                if d > 1e-10 and first is None:
                    first = (t, d)
                worst = max(worst, d)
            drift = cc.gas_mass_drift(w)
            last_ok = max([float(t) for t in pc.times(w1) if t in done and float(t) < float(first[0])]
                          or [0.0]) if first else None
            if first:
                print(f"    within 1e-10 up to t = {last_ok} s ({last_ok / cc.T_SCALE:.1f} sqrt(h0/g))")
            report(f"{n} processors", worst <= 1e-10 and drift <= 1e-12,
                   f"p, Ug, U within {worst:.1e} of their largest"
                   + (f" (first above 1e-10 at t = {first[0]} s: {first[1]:.1e})" if first else "")
                   + f"; gas mass {drift:.1e}; {secs[f'gas{n}']:.0f} s; {len(done)} of "
               f"{len(pc.times(w1))} writes")

    if "K" in only:
        print("\nK: cost per step, gas on against gas off (report; F.9 target 1.5x)")
        per = (secs["gas"] / nst["gas"]) / (secs["dry"] / nst["dry"])
        print(f"  {per:.2f}x ({secs['gas']:.0f} s for {nst['gas']} steps against {secs['dry']:.0f} s "
              f"for {nst['dry']} steps, {cc.N} cells, serial)")

    print(f"\ntest_fluidised_collapse: {'PASS' if all(res) else 'FAIL'} ({sum(res)}/{len(res)}), "
          f"{time.perf_counter() - t_all:.0f} s")
    return 0 if all(res) else 1


if __name__ == "__main__":
    only = "ABCK"
    if "--only" in sys.argv:
        only = sys.argv[sys.argv.index("--only") + 1].replace(",", "")
    sys.exit(main(sys.argv[1], only))
