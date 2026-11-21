"""WP K.1 done-when checks: the compaction cap at one point (decision 0113).
Seconds.

capCheck.C (g++ against the kernel): Mohr-Coulomb and Drucker-Prager with
capLambda lambda* = 0.01, E = 50 MPa, nu = 0.3, phi 30, from an isotropic
1 kPa, in isotropic compression and in an oedometer (only eps_yy), loaded to
a mean stress of 100 kPa and unloaded to 10 kPa.

[1] Loading follows the normal compression line eps_v = lambda* ln(P/P0)
    (P0 = 1 kPa, normally consolidated at the start) within 1 % of the
    strain, and p_c = P.
[2] Unloading is elastic: P falls by K d eps_v (K = E/(3(1 - 2 nu))) within
    1 %, and p_c stays at its largest value.

    python3 cap_checks.py WORKDIR
"""

import math
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
LAM, E, NU, P0 = 0.01, 5e7, 0.3, 1000.0
K = E / (3 * (1 - 2 * NU))
fails = []


def report(name, ok, detail):
    print(f"  {name}: {'PASS' if ok else 'FAIL'} ({detail})", flush=True)
    if not ok:
        fails.append(name)


def main():
    work = os.path.abspath(sys.argv[1])
    os.makedirs(work, exist_ok=True)
    k = os.path.join(ROOT, "src", "celConstitutive", "kernel")
    exe = os.path.join(work, "capCheck")
    srcs = [os.path.join(k, f) for f in sorted(os.listdir(k)) if f.endswith(".C")]
    subprocess.run(["g++", "-std=c++17", "-O2", "-Wall", "-Wextra", "-Werror", f"-I{k}",
                    os.path.join(HERE, "capCheck.C"), *srcs, "-o", exe], check=True)
    out = subprocess.run([exe], check=True, capture_output=True, text=True).stdout
    rows = {}
    for line in filter(None, out.split("\n")):
        m, t, ph, ev, P, pc = line.split()
        rows.setdefault((m, t), []).append((int(ph), float(ev), float(P), float(pc)))
    for (m, t), r in rows.items():
        a = np.array(r)
        load, unl = a[a[:, 0] == 0], a[a[:, 0] == 1]
        print(f"[{m} {t}] loaded to {load[-1, 2] / 1e3:.1f} kPa (eps_v {load[-1, 1]:.4f}), "
              f"unloaded to {unl[-1, 2] / 1e3:.1f} kPa", flush=True)
        ref = LAM * np.log(load[:, 2] / P0)
        err = np.abs(load[:, 1] - ref).max() / ref.max()
        pce = np.abs(load[:, 3] / load[:, 2] - 1).max()
        report("loading on lambda*", err <= 0.01 and pce <= 0.01 and len(load) > 100,
               f"largest strain error {err * 100:.3f} % of {ref.max():.4f}; p_c/P - 1 up to {pce:.1e}")
        dP = load[-1, 2] - unl[-1, 2]
        dE = load[-1, 1] - unl[-1, 1]
        kerr = abs(dP / (K * dE) - 1)
        pcs = abs(unl[-1, 3] / load[-1, 3] - 1)
        report("unloading elastic", kerr <= 0.01 and pcs <= 1e-12,
               f"dP/(K deps_v) - 1 = {kerr:.1e}; p_c kept ({pcs:.0e})")
    print(f"{len(fails)} failures")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
