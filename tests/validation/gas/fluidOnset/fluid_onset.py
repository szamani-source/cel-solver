"""Gas validation pack, test 2 (WP F.7, decisions 0125 item 9, 0132):
onset of fluidisation against Wen and Yu.

A bed of quartz sand (d 0.2 mm, grain density 2650 kg/m^3, packed at gas
fraction 0.38, 5 cm deep, 2 cm wide, 8 x 20 cells, 20 rows of void above)
rests on the floor of a column with smooth walls and an open top. Air is
blown in through the floor at a superficial velocity ramped from 0 to 1.4
times Wen and Yu's U_mf over 0.28 s (slow: the bed's pressure settles in
about 4 ms). Plain Mohr-Coulomb grains (phi 30, c 0), model gasFlow with the
pore gas momentum, drag ergunWenYu (Sia, 2026-11-24).

U_mf from the run: the air speed at which the pressure drop over the bed
(the floor's pressure, fitted u + u^2 below 0.85 of the weight) reaches the
bed's weight per area. Reference: Wen and Yu (1966),
Re_mf = sqrt(33.7^2 + 0.0408 Ar) - 33.7, tools/gas_reference.py. Their
constants stand for packings near e_mf = 0.38; U_mf of an Ergun-type law
moves about 9 % per 0.01 of e_mf.

Checks: U_mf within 10 % of Wen-Yu on 1, 2 and 4 processors (2 x 1: a
vertical cut through the bed and the free surface; 2 x 2: also a horizontal
cut at the free surface); 2 and 4 processors match serial to 1e-10 of each
field's largest value (p, Ug, U) up to the onset (after it the lifting bed
parts by round-off amplified to about 1e-2: reported, F.9); gas mass less
the inflow constant to 1e-12. Writes times.csv and results.md.

    python3 fluid_onset.py [WORKDIR]      # default ./run
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..", "..", "..")
sys.path.insert(0, os.path.join(ROOT, "tests", "gas"))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import fluid_cases as fl  # noqa: E402
import gas_flow_checks as gf  # noqa: E402
import pore_cases as pc  # noqa: E402

NX = 8
DECOMP = {1: (1, 1, 1), 2: (2, 1, 1), 4: (2, 2, 1)}
N = NX * (fl.NB + fl.NV)


def main():
    work = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "run")
    os.makedirs(work, exist_ok=True)
    rows, ok, serial = [], True, None
    for n, dec in DECOMP.items():
        path = f"{work}/np{n}"
        fl.make_bed(path, decomp=dec, nx=NX, extra_gas="pTolerance 1e-12;\n")
        secs = pc.run(path, n)
        u_dp, u_lift, a = fl.onset(path, nx=NX)
        drift = gf.drift(path)
        steps = gf.steps(path)
        if serial is None:
            serial = path
            match, after = 0.0, 0.0
        else:
            b, _, aft = fl.compare_runs(serial, path, N, u_dp / fl.U_END * fl.T_RAMP)
            match, after = max(b.values()), max(aft.values())
        err = u_dp / fl.U_WY - 1
        good = abs(err) <= 0.10 and match <= 1e-10 and drift <= 1e-12
        ok = ok and good
        rows.append((n, steps, secs, u_dp, fl.U_WY, err, u_lift, match, after, drift, good))
        print(f"np {n}: U_mf {u_dp:.5f} m/s, Wen-Yu {fl.U_WY:.5f} m/s, {100 * err:+.2f} %; lift from "
              f"{u_lift:.5f} m/s; vs serial {match:.0e} to the onset ({after:.0e} after); mass {drift:.0e}; {steps} steps {secs:.1f} s: "
              f"{'PASS' if good else 'FAIL'}", flush=True)

    with open(os.path.join(HERE, "times.csv"), "w") as f:
        f.write("nprocs,steps,wall_s,U_mf_m_per_s,wen_yu_m_per_s,error,U_lift_m_per_s,"
                "vs_serial_to_onset,vs_serial_after,mass_drift,pass\n")
        for r in rows:
            f.write("%d,%d,%.3f,%.6f,%.6f,%.4e,%.6f,%.1e,%.1e,%.1e,%s\n" % r)
    with open(os.path.join(HERE, "results.md"), "w") as f:
        f.write("# Gas validation pack test 2: onset of fluidisation against Wen and Yu\n\n"
                f"Quartz sand d {fl.D_BED * 1e3:g} mm, packed at gas fraction {fl.E_MF}, bed {fl.HB} m "
                f"deep ({NX} x {fl.NB} cells, {fl.NV} void rows above), air at {pc.T} K ramped to "
                f"{fl.U_END:.4f} m/s over {fl.T_RAMP} s; model gasFlow with the pore gas momentum, "
                f"drag ergunWenYu. Drag law's own U_mf at this packing {fl.U_MODEL:.5f} m/s. "
                "Written by `fluid_onset.py`.\n\n"
                "| processors | steps | wall [s] | U_mf [m/s] | Wen-Yu [m/s] | error | lift from [m/s] | "
                "vs serial to the onset | after (report) | gas mass |\n"
                "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |\n")
        for r in rows:
            f.write("| %d | %d | %.1f | %.5f | %.5f | %+.2f %% | %.5f | %.0e | %.0e | %.0e |\n"
                    % (r[0], r[1], r[2], r[3], r[4], 100 * r[5], r[6], r[7], r[8], r[9]))
        f.write(f"\nAll checks {'pass' if ok else 'FAIL'}: within 10 % of Wen and Yu, processor "
                "counts agree to 1e-10 of each field's largest value (p, Ug, U) up to the onset, "
                "gas mass to 1e-12. After the onset the lifting bed parts serial from parallel "
                "(round-off amplified; reported, F.9, decision 0132).\n")
    print(f"\nFluid onset: {'ALL PASS' if ok else 'FAILED'}; times.csv and results.md written")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
