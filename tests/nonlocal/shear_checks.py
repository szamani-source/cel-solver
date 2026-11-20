"""WP N.2 done-when checks: a granular layer sheared by a moving lid under
gravity (shear_layer.py, decision 0111). About 3 minutes.

[1] Nonlocal law on 20, 40 and 80 cells (2, 1 and 0.5 d): the band width
    (height over which u goes from 0.1 to 0.9 V) on 1 and 0.5 d within 5 %
    of each other and of the 1D nonlocal solution for the run's own P(y)
    (tests/nonlocal/shear_layer.py reference), and the velocity profile
    within 5 % of V of it. The law's length is a few d here (the band is
    about 3 d), so 2 d cells are too coarse to resolve it: they are run and
    reported, not checked.
[2] Local law (no nonlocalA) on 20 and 40 cells, for contrast: the lid
    slips over the top cell and the band follows the cell. Reported only.

    python3 shear_checks.py WORKDIR
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import shear_layer as sl  # noqa: E402

fails = []


def report(name, ok, detail):
    print(f"  {name}: {'PASS' if ok else 'FAIL'} ({detail})", flush=True)
    if not ok:
        fails.append(name)


def main():
    work = os.path.abspath(sys.argv[1])
    os.makedirs(work, exist_ok=True)
    print("[1] Nonlocal law, 2 / 1 / 0.5 d cells", flush=True)
    sys.argv = [sys.argv[0], work, "nonlocal", "20", "40", "80"]
    r = sl.main()
    d = sl.LAW["d"]
    b40, ref40, e40 = r[40]
    b80, ref80, e80 = r[80]
    report("band converges", abs(b80 / b40 - 1) <= 0.05,
           f"{r[20][0] / d:.2f} / {b40 / d:.2f} / {b80 / d:.2f} d (1 and 0.5 d within 5 %)")
    report("band as 1D nonlocal", max(abs(b40 / ref40 - 1), abs(b80 / ref80 - 1)) <= 0.05,
           f"{b40 / d:.2f} against {ref40 / d:.2f} d, {b80 / d:.2f} against {ref80 / d:.2f} d")
    report("profile as 1D nonlocal", max(e40, e80) <= 0.05,
           f"largest difference {e40 * 100:.1f} / {e80 * 100:.1f} % of V (<= 5 %)")
    print("[2] Local law, 2 / 1 d cells (reported)", flush=True)
    sys.argv = [sys.argv[0], work, "local", "20", "40"]
    sl.main()
    print(f"{len(fails)} failures")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
