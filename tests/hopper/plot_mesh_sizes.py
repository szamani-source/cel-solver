"""Cell-size picture of a Couto case mesh (decision 0115 addendum):
    python3 tests/hopper/plot_mesh_sizes.py CASE PNG TITLE   (from the repo root)
"""
import math, sys, numpy as np
sys.path.insert(0, "tests/angled"); sys.path.insert(0, "tests/hopper")
import angled_checks as ac, hopper_case as hc, buffered_case as bc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
p, png, title = sys.argv[1], sys.argv[2], sys.argv[3]
m = f"{p}/constant/polyMesh"
pts = ac.foam_list(f"{m}/points", "point"); faces = ac.foam_list(f"{m}/faces", "face")
t = math.radians(hc.THETA)
polys = [np.column_stack([pts[f][:, 0] / math.cos(t), pts[f][:, 1]])
         for f in faces if np.all(pts[f][:, 2] <= 1e-12) and len(f) == 4]
def size(q):
    a = q[[1, 2, 3, 0]] - q
    L = np.linalg.norm(a, axis=1)
    return 1000 * min(L.max(), max(L[0], L[2]), max(L[1], L[3]))
val = np.array([1000 * min(max(np.linalg.norm(q[1]-q[0]), np.linalg.norm(q[3]-q[2])),
                           max(np.linalg.norm(q[2]-q[1]), np.linalg.norm(q[0]-q[3]))) for q in polys])
c = dict(hc.COUTO)
poly = bc.shell_polygon(c, y_top=3.0)
fig, axs = plt.subplots(1, 5, figsize=(22, 9), gridspec_kw=dict(width_ratios=[1.1, 1, 1, 1, 1]))
views = [((-0.02, 0.6), (-0.05, 2.85), "Whole silo"),
         ((0.10, 0.24), (-0.03, 0.08), "Outlet edge, close up"),
         ((0.16, 0.36), (-0.03, 0.14), "Outlet edge and hopper wall"),
         ((0.36, 0.56), (0.40, 0.62), "Transition (wall)"),
         ((0.0, 0.56), (1.40, 1.66), "Upper cylinder (axis to wall)")]
for ax, (xl, yl, ttl) in zip(axs, views):
    pc = PolyCollection(polys, array=val, cmap="Spectral", clim=(1, 16),
                        edgecolors="#ffffff", linewidths=0.35)
    ax.add_collection(pc)
    pp = np.vstack([poly, poly[:1]])
    ax.plot(pp[:, 0], pp[:, 1], color="#000000", lw=1.0, label="STL shell (20 mm)")
    ax.set_xlim(*xl); ax.set_ylim(*yl); ax.set_aspect("equal"); ax.set_title(ttl); ax.set_xlabel("r (m)")
axs[0].set_ylabel("y (m)"); axs[0].legend(loc="upper left", fontsize=8)
fig.colorbar(pc, ax=axs, shrink=0.6, label="cell size, smaller side (mm)")
fig.suptitle(title)
fig.savefig(png, dpi=130); print(png, len(polys), np.percentile(val,[0,50,100]))
