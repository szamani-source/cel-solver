"""WP 3.6f: scatter of the T0 discharge rate (decision 0046); run by
tests/manual/test_discharge_scatter.sh.

A plastic flow parts by round-off growth (decision 0032), so one run gives
one sample of the rate. This runs each E several times with a slightly
different time step (maxCo 0.40, 0.38, 0.36, ...: the same physics, a
different round-off history) and prints each rate, the mean and the range.
Material as tests/discharge/discharge_checks.py --mui.

    python3 scatter.py WORKDIR [--cell 0.01] [--runs 3] [--jobs 3] [E_MPa ...]

Extra material entries come from CEL_EXTRA_MATERIAL, as in every driver.
"""

import os
import sys
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import discharge_checks as dc  # noqa: E402

oc = dc.oc


def one(job):
    work, E_mpa, co, cell = job
    oc.CELL = oc.vc.CELL = cell
    oc.MAXCO = co
    sub = f"{work}/co{co:g}"
    os.makedirs(sub, exist_ok=True)
    r = dc.case(sub, E_mpa)
    return E_mpa, co, r["rate"], r["C"], r["sec"]


def main():
    args = sys.argv[1:]
    work = args.pop(0)
    opt = {"--cell": 0.01, "--runs": 3, "--jobs": 3}
    for k in list(opt):
        if k in args:
            i = args.index(k)
            opt[k] = type(opt[k])(args[i + 1])
            del args[i:i + 2]
    elist = [float(a) for a in args] or dc.E_LIST
    os.environ["CEL_EXTRA_MATERIAL"] = dc.MUI + os.environ.get("CEL_EXTRA_MATERIAL", "")
    cos = [round(0.40 - 0.02 * k, 2) for k in range(opt["--runs"])]
    jobs = [(f"{work}/E{e:g}", e, c, opt["--cell"]) for e in elist for c in cos]
    with Pool(opt["--jobs"]) as p:
        res = p.map(one, jobs)
    print(f"T0 tank, mu(I) + flow impedance, {1000 * opt['--cell']:g} mm cells, "
          f"maxCo {', '.join(f'{c:g}' for c in cos)}"
          + (f"; extra: {os.environ.get('CEL_EXTRA_MATERIAL', '')[len(dc.MUI):].strip()}"
             if os.environ["CEL_EXTRA_MATERIAL"] != dc.MUI else ""))
    for e in elist:
        r = np.array([x[2] for x in res if x[0] == e])
        print(f"  E = {e:g} MPa: " + " ".join(f"{v:.1f}" for v in r)
              + f" kg/s/m; mean {r.mean():.1f}, range {r.min():.1f} to {r.max():.1f} "
              f"({100 * (r.max() - r.min()) / r.mean():.0f} % of the mean)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
