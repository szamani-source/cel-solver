"""Gas validation pack, test 3 (WP F.3, decisions 0125 item 9, 0128):
pore-pressure decay in a suddenly loaded fine-powder column against the
exact diffusion (Terzaghi) solution.

Fine quartz powder: d 0.1 mm, grain density 2650 kg/m^3, gas fraction
0.40, linear elastic skeleton E 1 MPa, nu 0.3 (constrained modulus M
1.346 MPa); column 1 m tall, 80 cells of 12.5 mm, one wide, smooth sealed
walls and base, top drained (p = 0 gauge) and loaded with 2 kPa from
t = 0; no gravity; air at 293.15 K and 101325 Pa; drag ergunWenYu;
velocity damping 50 1/s against the elastic waves of the sudden load.
The solid and the gas are coupled both ways (effective stress).

Reference: Terzaghi's one-dimensional consolidation with the gas storage,
c = k/(mu (e/p0 + 1/M)) = 0.139 m^2/s (k Ergun's permeability), starting
from the undrained pressure p_i = q (1/M)/(e/p0 + 1/M) = 316.7 Pa:
degree U(T) = 1 - mean(p)/p_i and the profile p(z, T)/p_i (series,
tools/gas_reference.py), T = c t / H^2, run to T = 1 (7.2 s).

Checks: U within 3 % of Terzaghi at T = 0.05, 0.1, 0.2, 0.5 and 1; the
profile within 0.03 p_i at T >= 0.1; 2 and 4 processors (cuts across the
column) match serial to 1e-10 of p_i; gas mass less the
boundary flow constant to 1e-12. Writes times.csv and results.md.

A column 4 cells wide cut also along its length (2 x 2) parts from serial by
1e-9 of p_i: the solid itself does (4e-6 of the load with gas off), a
known serial/parallel issue of the solid (decision 0128, STATUS known
issues), not of the gas.

    python3 terzaghi_column.py [WORKDIR]      # default ./run
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..", "..", "..")
sys.path.insert(0, os.path.join(ROOT, "tests", "gas"))
import consolidation_cases as cc  # noqa: E402
import pore_cases as pc  # noqa: E402

NX, NY = 1, 80
TS = (0.05, 0.1, 0.2, 0.5, 1.0)
DECOMP = {1: (1, 1, 1), 2: (1, 2, 1), 4: (1, 4, 1)}


def main():
    work = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "run")
    os.makedirs(work, exist_ok=True)
    rows, table, ok, serial = [], [], True, None
    for n, dec in DECOMP.items():
        path = f"{work}/np{n}"
        cc.make(path, nx=NX, ny=NY, end=cc.t_of(1.0), write_interval=cc.t_of(0.05), decomp=dec)
        secs = pc.run(path, n)
        ts = {}
        for t in pc.times(path):
            T = float(t) * cc.C_V / cc.H ** 2
            for x in TS:
                if abs(T - x) < 1e-6:
                    ts[x] = t
        M, N, ns, _ = pc.mass_log(path)
        drift = float(np.abs(M - N - (M[0] - N[0])).max() / M[0])
        fields = {T: pc.field(path, ts[T], "p", NX * NY) for T in TS}
        if serial is None:
            serial = fields
        match = max(float(np.abs(fields[T] - serial[T]).max()) for T in TS) / cc.P_I
        good = match <= 1e-10 and drift <= 1e-12
        worst = 0.0
        for T in TS:
            p = fields[T].reshape(NY, NX).mean(axis=1) / cc.P_I
            U = 1 - float(p.mean())
            Ut = cc.theory(T)
            perr = float(np.abs(p - cc.profile_theory(NY, T)).max())
            g = abs(U - Ut) <= 0.03 * Ut and (T < 0.1 or perr <= 0.03)
            good = good and g
            worst = max(worst, abs(U - Ut) / Ut)
            if n == 1:
                table.append((T, U, Ut, (U - Ut) / Ut, perr))
        ok = ok and good
        rows.append((n, len(M), secs, worst, match, drift, ns.mean(), good))
        print(f"np {n}: U within {100 * worst:.2f} % of Terzaghi; vs serial {match:.0e}; "
              f"mass {drift:.0e}; {len(M)} steps {secs:.1f} s: {'PASS' if good else 'FAIL'}",
              flush=True)

    with open(os.path.join(HERE, "times.csv"), "w") as f:
        f.write("nprocs,steps,wall_s,worst_U_error,vs_serial,mass_drift,solves_per_step,pass\n")
        for r in rows:
            f.write("%d,%d,%.3f,%.3e,%.1e,%.1e,%.2f,%s\n" % r)
    with open(os.path.join(HERE, "results.md"), "w") as f:
        f.write("# Gas validation pack test 3: pore-pressure decay against Terzaghi\n\n"
                f"Fine quartz powder d {cc.D_GRAIN * 1e3:g} mm, gas fraction {cc.E_GAS}, "
                f"E {cc.E_MOD / 1e6:g} MPa, nu {cc.NU} (M {cc.M_MOD / 1e6:.4g} MPa); column "
                f"{cc.H} m ({NX} x {NY} cells), drained top loaded with {cc.Q:g} Pa at t = 0; "
                f"air at {pc.T} K; drag ergunWenYu; damping {cc.DAMPING:g} 1/s. "
                f"c = {cc.C_V:.4f} m^2/s, p_i = {cc.P_I:.2f} Pa. Written by "
                "`terzaghi_column.py`.\n\n"
                "Serial run against the series:\n\n"
                "| T | U | Terzaghi | error | profile, largest error / p_i |\n"
                "| --- | --- | --- | --- | --- |\n")
        for T, U, Ut, e, perr in table:
            f.write(f"| {T:g} | {U:.4f} | {Ut:.4f} | {100 * e:+.2f} % | {perr:.3f} |\n")
        f.write("\n| processors | steps | wall [s] | worst U error | vs serial | gas mass | "
                "solves per step |\n| --- | --- | --- | --- | --- | --- | --- |\n")
        for r in rows:
            f.write("| %d | %d | %.1f | %.2f %% | %.0e | %.0e | %.2f |\n"
                    % (r[0], r[1], r[2], 100 * r[3], r[4], r[5], r[6]))
        f.write(f"\nAll checks {'pass' if ok else 'FAIL'}: U within 3 % of Terzaghi, profiles "
                "within 0.03 p_i from T = 0.1, processor counts agree to 1e-10 of p_i, gas mass "
                "to 1e-12. Before T = 0.05 the elastic waves of the sudden load (not in the "
                "quasi-static theory) are still being damped.\n")
    print(f"\nTerzaghi column: {'ALL PASS' if ok else 'FAILED'}; times.csv and results.md written")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
