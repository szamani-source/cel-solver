"""Gas validation pack, test 1 (WP F.2, decisions 0125 item 9, 0127):
packed-bed pressure drop against Ergun.

Air is blown up through a fixed column of quartz sand: d 0.5 mm, grain
density 2650 kg/m^3, gas fraction 0.40; column 0.05 m wide, 0.3 m tall,
10 x 60 cells (2D), sides sealed, top open (p = 0 gauge), inflow at the
bottom at 0.02, 0.05, 0.1 and 0.2 m/s (superficial, at 101325 Pa and
293.15 K). The sand is frozen (no gravity), so this is the rigid-bed
pressure drop; minimum fluidisation (Wen-Yu, about 0.2 m/s for this sand)
is pack test 2 (F.7). Drag law ergunWenYu. Each run goes to 0.1 s (about
five times the pressure relaxation time L^2/c = 0.02 s).

Reference: Ergun (1952) integrated up the column for an isothermal ideal
gas, p_in^2 = p_out^2 + 2 R T (a m + b m^2) L, with a = 150 mu (1-e)^2 /
(e^3 d^2), b = 1.75 (1-e)/(e^3 d), m = rho_g U the mass flux.

Checks: the inlet pressure within 1 % of Ergun at every speed; 2 and 4
processors (cuts through the bed) match serial to 1e-10 of the pressure
drop; gas mass balance (mass - inflow) constant to 1e-12. Writes
times.csv (one row per run) and results.md.

    python3 ergun_bed.py [WORKDIR]      # default ./run
"""

import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..", "..", "..")
sys.path.insert(0, os.path.join(ROOT, "tests", "gas"))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import gas_reference as gr  # noqa: E402
import pore_cases as pc  # noqa: E402

D, E_BED, W, L, NX, NY = 0.5e-3, 0.40, 0.05, 0.3, 10, 60
SPEEDS = (0.02, 0.05, 0.1, 0.2)
DECOMP = {1: (1, 1, 1), 2: (1, 2, 1), 4: (2, 2, 1)}
MU, RG = gr.air_viscosity(pc.T), gr.air_density(pc.P0, pc.T)
RT = gr.R_AIR * pc.T


def ergun_inlet(U):
    m = U * RG
    a = 150.0 * MU * (1 - E_BED) ** 2 / (E_BED ** 3 * D * D)
    b = 1.75 * (1 - E_BED) / (E_BED ** 3 * D)
    return math.sqrt(pc.P0 ** 2 + 2 * RT * (a * m + b * m * m) * L) - pc.P0


def main():
    work = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "run")
    os.makedirs(work, exist_ok=True)
    rows, ok = [], True
    for U in SPEEDS:
        ref = ergun_inlet(U)
        serial = None
        for n, dec in DECOMP.items():
            path = f"{work}/U{U}_np{n}"
            pc.make_case(path, size=(W, L), cells=(NX, NY), e=E_BED, d=D, end=0.1,
                         write_interval=0.05,
                         p_patches={"top": "type fixedValue; value uniform 0;",
                                    "bottom": "type fixedGradient; gradient uniform 0;"},
                         inflow={"bottom": U}, decomp=dec)
            secs = pc.run(path, n)
            t = pc.times(path)[-1]
            p = pc.field(path, t, "p", NX * NY)
            g = pc.patch_gradient(path, t, "p", "bottom")
            pin = float(np.mean(p[:NX] + np.atleast_1d(g) * (L / NY) / 2))
            M, N, ns, _ = pc.mass_log(path)
            drift = np.abs(M - N - (M[0] - N[0])).max() / M[0]
            steps = len(M)
            if serial is None:
                serial = p
            match = np.abs(p - serial).max() / ref
            err = pin / ref - 1
            good = abs(err) <= 0.01 and match <= 1e-10 and drift <= 1e-12
            ok = ok and good
            rows.append((U, n, steps, secs, pin, ref, err, match, drift, ns.mean(), good))
            print(f"U {U:5.2f} m/s  np {n}: inlet {pin:9.3f} Pa, Ergun {ref:9.3f} Pa, "
                  f"{100 * err:+.4f} %; vs serial {match:.0e}; mass {drift:.0e}; "
                  f"{steps} steps {secs:.1f} s: {'PASS' if good else 'FAIL'}", flush=True)

    with open(os.path.join(HERE, "times.csv"), "w") as f:
        f.write("U_m_per_s,nprocs,steps,wall_s,dp_Pa,ergun_Pa,error,vs_serial,mass_drift,"
                "solves_per_step,pass\n")
        for r in rows:
            f.write("%g,%d,%d,%.3f,%.6f,%.6f,%.3e,%.1e,%.1e,%.2f,%s\n" % r)
    with open(os.path.join(HERE, "results.md"), "w") as f:
        f.write("# Gas validation pack test 1: packed-bed pressure drop against Ergun\n\n"
                f"Quartz sand d {D * 1e3:g} mm, gas fraction {E_BED}, column {L} m "
                f"({NX} x {NY} cells), air at {pc.T} K; drag ergunWenYu. "
                "Written by `ergun_bed.py`.\n\n"
                "| U [m/s] | processors | steps | wall [s] | inlet [Pa] | Ergun [Pa] | error | "
                "vs serial | gas mass |\n| --- | --- | --- | --- | --- | --- | --- | --- | --- |\n")
        for r in rows:
            f.write("| %g | %d | %d | %.1f | %.3f | %.3f | %+.4f %% | %.0e | %.0e |\n"
                    % (r[0], r[1], r[2], r[3], r[4], r[5], 100 * r[6], r[7], r[8]))
        f.write(f"\nAll checks {'pass' if ok else 'FAIL'}: within 1 % of Ergun, processor "
                "counts agree to 1e-10 of the pressure drop, gas mass to 1e-12.\n")
    print(f"\nErgun bed: {'ALL PASS' if ok else 'FAILED'}; times.csv and results.md written")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
