"""WP 3.6s: are the shear bands of the column collapse mesh dependent?

Collapse a = 0.5 at E = 1 MPa (the tests/test_collapse.sh set-up, shear cap
and rule 9 on, decision 0032), with the plain T0 Mohr-Coulomb and with mu(I)
(decision 0037 placeholders), on 10, 5 and 2.5 mm cells (or the sizes given).

Measured:
- band width: at t = 0.1 and 0.2 s (failure wedge sliding), the equivalent
  shear rate gdot = sqrt((ux,x - uy,y)^2 + (ux,y + uy,x)^2) on full cells;
  in each row between 20 and 60 % of H_i, the run of cells around the row's
  peak with gdot >= half the peak (FWHM, horizontal), averaged over the rows.
  A band that is a fixed number of cells wide (width in m halves with the
  cell) is mesh dependent; a band of fixed width in m is converged.
- run-out (L_f - L_i)/L_i and final height at t = 0.5 and 1 s, as in
  collapse_checks.measure (front: half a cell of mass, so it depends on the
  cell size), and with the front where the layer is 5 mm thick.

    python3 shear_band.py WORKDIR [cell_mm ...]

Runs the cases on up to 4 processes at once (one case per process).
"""

import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "void"))
import collapse_checks as cc  # noqa: E402
import void_checks as vc  # noqa: E402

A, E_MPA, T_MAX = 0.5, 1.0, 1.0
SNAP = [0.1, 0.2]
MU_I = "mu2 0.961813;\nI0 0.279;\nd 0.001;\nrhoS 2650;\n"   # decision 0037
# WP N.3 (decision 0111): nonlocal fluidity on the mu(I) set, with the grain
# size d set so the band (a few d) is resolved on 5 and 2.5 mm cells; A of
# Henann and Kamrin (2013), glass beads
NL_D = float(os.environ.get("CEL_NL_D", "0.005"))
NONLOCAL = (MU_I.replace("d 0.001;", f"d {NL_D!r};")
            + "nonlocalA 0.48;\nflowImpedance yes;\n")
LAWS = os.environ.get("CEL_LAWS", "MC muI").split()


def set_cell(h):
    vc.CELL = h
    cc.CELL = h


def band_width(cells, nx, ny, h, H_i):
    c = np.loadtxt(cells)
    i = np.floor(c[:, 0] / h + 1e-9).astype(int)
    j = np.floor(c[:, 1] / h + 1e-9).astype(int)
    al = np.zeros((ny, nx))
    ux = np.zeros((ny, nx))
    uy = np.zeros((ny, nx))
    al[j, i], ux[j, i], uy[j, i] = c[:, 2], c[:, 4], c[:, 5]
    full = al > 0.95
    ok = full.copy()
    ok[1:-1, 1:-1] &= full[2:, 1:-1] & full[:-2, 1:-1] & full[1:-1, 2:] & full[1:-1, :-2]
    ok[0, :] = ok[-1, :] = False
    ok[:, 0] = ok[:, -1] = False
    g = np.zeros_like(ux)
    dxx = (ux[1:-1, 2:] - ux[1:-1, :-2]) / (2 * h)
    dyx = (ux[2:, 1:-1] - ux[:-2, 1:-1]) / (2 * h)
    dxy = (uy[1:-1, 2:] - uy[1:-1, :-2]) / (2 * h)
    dyy = (uy[2:, 1:-1] - uy[:-2, 1:-1]) / (2 * h)
    g[1:-1, 1:-1] = np.sqrt((dxx - dyy) ** 2 + (dyx + dxy) ** 2)
    g[~ok] = 0.0
    widths, peaks = [], []
    for jj in range(ny):
        y = (jj + 0.5) * h
        if not (0.2 * H_i <= y <= 0.6 * H_i):
            continue
        row = g[jj]
        k = int(row.argmax())
        if row[k] <= 0.0:
            continue
        half = 0.5 * row[k]
        lo = k
        while lo > 0 and row[lo - 1] >= half:
            lo -= 1
        hi = k
        while hi < nx - 1 and row[hi + 1] >= half:
            hi += 1
        widths.append((hi - lo + 1) * h)
        peaks.append(row[k])
    return (float(np.mean(widths)) if widths else float("nan"),
            float(np.mean(peaks)) if peaks else float("nan"))


def one(args):
    work, h, law = args
    set_cell(h)
    L_i, H_i, nx, ny = cc.geometry(A)
    name = f"band_{law}_{h * 1000:g}mm"
    path = cc.setup(work, A, E_MPA, T_MAX, name)
    dct = f"{path}/system/voidCheckDict"
    txt = re.sub(r"writeTimes[^;]*;", "writeTimes (" + " ".join(str(t) for t in SNAP + [0.5])
                 + ");", open(dct).read())
    open(dct, "w").write(txt)
    if law in ("muI", "NL"):
        mp = f"{path}/constant/materialProperties"
        add = MU_I if law == "muI" else NONLOCAL
        txt = open(mp).read().replace("sigmaT 0;\n", "sigmaT 0;\n" + add)
        open(mp, "w").write(txt)
    sec = vc.run(path)
    out = dict(law=law, h=h, nx=nx, ny=ny, seconds=sec)
    for t in SNAP:
        out[f"w{t}"], out[f"g{t}"] = band_width(f"{path}/cells_{t:g}.dat", nx, ny, h, H_i)
    r05 = cc.measure(f"{path}/cells_0.5.dat", nx, ny, L_i)
    r1 = cc.measure(f"{path}/cells.dat", nx, ny, L_i)
    out.update(run05=r05["runout"], run1=r1["runout"], H_f=r1["H_f"], slope=r1["slope"])
    for key, r in (("f05", r05), ("f1", r1)):     # front where the layer is 5 mm thick
        k = np.where(r["hc"] >= 0.005)[0]
        out[key] = ((k.max() + 1) * h - L_i) / L_i
    hist = np.loadtxt(f"{path}/history.dat")
    out["dm"] = abs(hist[-1, cc.MASS] / hist[0, cc.MASS] - 1.0)
    return out


def main():
    work = sys.argv[1]
    sizes = [float(s) / 1000 for s in sys.argv[2:]] or [0.01, 0.005, 0.0025]
    jobs = [(work, h, law) for h in sizes for law in LAWS]
    t0 = time.perf_counter()
    with ProcessPoolExecutor(max_workers=min(4, len(jobs))) as ex:
        res = list(ex.map(one, jobs))
    print(f"Column collapse a = {A}, E = {E_MPA} MPa, shear cap and rule 9 on, t_max {T_MAX} s")
    print("law  cell   band width at 0.1 s / 0.2 s     peak gdot 0.1/0.2  run-out 0.5 s / 1 s"
          "  (5 mm front)    H_f     slope  mass    wall")
    for r in res:
        hmm = r["h"] * 1000
        print(f"{r['law']:4s} {hmm:4.1f}mm "
              f"{r['w0.1'] * 1000:5.1f} / {r['w0.2'] * 1000:5.1f} mm "
              f"({r['w0.1'] / r['h']:4.1f} / {r['w0.2'] / r['h']:4.1f} cells)  "
              f"{r['g0.1']:6.1f} / {r['g0.2']:6.1f} 1/s  "
              f"{r['run05']:.3f} / {r['run1']:.3f}  ({r['f05']:.3f} / {r['f1']:.3f})  {r['H_f']:.4f}  {r['slope']:4.1f}  "
              f"{r['dm']:.0e} {r['seconds']:5.0f} s")
    print(f"total {time.perf_counter() - t0:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
