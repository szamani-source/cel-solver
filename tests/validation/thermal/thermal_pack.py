"""Thermal validation pack (WP T.8, decision 0136).

Three cases from the done-when tests, each run on 1, 2 and 4 processors:

schumann  hot air into a cold packed column (T.4): theta_g and theta_s
          against Schumann (1929) with Gunn's h; limit 1 % of the step.
graetz    laminar air between plates at a fixed wall temperature, in the
          void (T.5): developed Nu_Dh against 7.5407; limit 3 %.
movingBed ceramic particles between plates heated by a uniform flux, the
          case of Adapa et al. (2024) (T.6): mean h against their slug-flow
          solution with the near-wall gap (eq. 9); limit 5 %.

Every processor count must match serial to 1e-10 of the temperature
range. Writes times.csv (one row per run) and results.md.

    python3 thermal_pack.py [WORKDIR]      # default ./run
"""

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..", "..")
sys.path.insert(0, os.path.join(ROOT, "tests", "thermal"))
import conduction_checks as cc  # noqa: E402
import gas_heat_checks as gh  # noqa: E402
import moving_bed_checks as mb  # noqa: E402
import pore_cases as pc  # noqa: E402
import void_heat_checks as vh  # noqa: E402

DECOMP = {1: (1, 1, 1), 2: (2, 1, 1), 4: (4, 1, 1)}


def schumann(path, dec):
    gh.column(path, "twoTemperature", gh.T0 + gh.DT0, decomp=dec)
    g, s = gh.schumann_errors(path)
    return max(g, s), 0.01, f"theta_g {g * 100:.2f} %, theta_s {s * 100:.2f} %", \
        ("T", "Tg"), gh.END, gh.DT0


def graetz(path, dec):
    vh.make(path, dec)
    pc.run(path, dec[0] * dec[1] * dec[2])
    nu, _, _ = vh.nusselt(path)
    want = vh.tr.graetz_plates_nu()
    return abs(nu / want - 1), 0.03, f"Nu {nu:.4f} against {want:.4f}", ("Tg",), vh.END, vh.DT0


def moving_bed(path, dec):
    mb.make(path, dec)
    cc.run(path, dec[0] * dec[1] * dec[2])
    worst, h, hw = mb.metrics(path)
    rng = np.ptp(cc.read_internal(f"{path}/{mb.END:g}/T"))
    return abs(h / hw - 1), 0.05, f"h {h:.1f} against {hw:.1f} W/(m^2 K), local Nu {worst * 100:.2f} %", \
        ("T",), mb.END, rng


CASES = {"schumann": schumann, "graetz": graetz, "movingBed": moving_bed}


def main():
    work = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "run")
    os.makedirs(work, exist_ok=True)
    rows, ok = [], True
    for name, fn in CASES.items():
        serial = None
        for n, dec in DECOMP.items():
            p = os.path.join(work, f"{name}{n}")
            t0 = time.perf_counter()
            err, lim, detail, fields, t_end, scale = fn(p, dec)
            secs = time.perf_counter() - t0
            if serial is None:
                serial, vs = p, 0.0
            else:
                vs = max(np.abs(cc.read_internal(f"{p}/{t_end:g}/{f}")
                                - cc.read_internal(f"{serial}/{t_end:g}/{f}")).max()
                         for f in fields) / scale
            good = err < lim and vs < 1e-10
            ok &= good
            rows.append((name, n, secs, err, lim, vs, detail, good))
            print(f"{name} on {n}: error {err * 100:.3f} % (limit {lim * 100:g} %), "
                  f"vs serial {vs:.1e}, {secs:.0f} s: {'PASS' if good else 'FAIL'}", flush=True)
    with open(os.path.join(HERE, "times.csv"), "w") as f:
        f.write("case,nprocs,wall_s,error,limit,vs_serial,pass\n")
        for r in rows:
            f.write(f"{r[0]},{r[1]},{r[2]:.1f},{r[3]:.4e},{r[4]:g},{r[5]:.1e},{r[7]}\n")
    with open(os.path.join(HERE, "results.md"), "w") as f:
        f.write("# Thermal validation pack\n\nWritten by `thermal_pack.py` (decision 0136).\n\n"
                "| case | processors | wall [s] | result | error | limit | vs serial |\n"
                "| --- | --- | --- | --- | --- | --- | --- |\n")
        for r in rows:
            f.write(f"| {r[0]} | {r[1]} | {r[2]:.0f} | {r[6]} | {r[3] * 100:.2f} % | "
                    f"{r[4] * 100:g} % | {r[5]:.0e} |\n")
        f.write("\n" + ("All checks pass." if ok else "Some checks FAIL.") + "\n")
    print("thermal pack: " + ("PASS" if ok else "FAIL"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
