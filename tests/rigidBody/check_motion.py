"""WP 4.1 checks: rigid body poses against the exact prescribed paths.

Usage: check_motion.py <serial case> <parallel case>
"""
import math
import os
import sys

import numpy as np

serial, parallel = sys.argv[1], sys.argv[2]
fails = []


def check(name, err, tol):
    ok = err <= tol
    print(f"{'ok  ' if ok else 'FAIL'} {name}: {err:.3e} (limit {tol:g})")
    if not ok:
        fails.append(name)


def read_stl(path):
    """Vertices of an ASCII STL, rounded to float32 as OpenFOAM reads them."""
    v = [list(map(float, l.split()[1:4]))
         for l in open(path) if l.strip().startswith("vertex")]
    return np.array(v, dtype=np.float32).astype(float)


def read_dat(case, body):
    path = os.path.join(case, "postProcessing", "rigidBodies", body + ".dat")
    return np.loadtxt(path, comments="#", ndmin=2)


def Rz(a):
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


tri = os.path.join(serial, "constant", "triSurface")
cube0 = read_stl(os.path.join(tri, "cube.stl"))
wall0 = read_stl(os.path.join(tri, "hopperWall.stl"))
cube_c0 = cube0.mean(axis=0)   # centroid of a cube surface = vertex mean


# Exact paths: centre c(t), rotation R(t), velocity v(t), omega w(t)
def pusher(t):
    T, V = 0.005, 1.0
    dz = -V*t*t/(2*T) if t <= T else -V*(T/2 + (t - T))
    v = -V*min(t/T, 1.0)
    return (cube_c0 + [0, 0, dz], np.eye(3), np.array([0, 0, v]),
            np.zeros(3))


def spinner(t):
    c0 = np.array([0.4, 0.4, 0.7])
    return (c0 + [0.5*t, 0, 0], Rz(100*t), np.array([0.5, 0, 0]),
            np.array([0, 0, 100.0]))


def wall(t):
    return None, np.eye(3), np.zeros(3), np.zeros(3)


bodies = {"pusher": (pusher, cube0, cube_c0),
          "spinner": (spinner, cube0, np.array([0.4, 0.4, 0.7])),
          "hopperWall": (wall, wall0, None)}

for body, (path, X0, c0) in bodies.items():
    d = read_dat(serial, body)
    if len(d) < 3:
        fails.append(body + " write times")
    ec = eR = ev = 0.0
    for row in d:
        t = row[0]
        c, R, v, w = path(t)
        if c is not None:
            ec = max(ec, np.abs(row[1:4] - c).max())
        eR = max(eR, np.abs(row[4:13].reshape(3, 3) - R).max())
        ev = max(ev, np.abs(row[13:16] - v).max(), np.abs(row[16:19] - w).max())
    print(f"{body}: {len(d)} write times, last t = {d[-1, 0]:.6g}")
    if c0 is not None:
        check(f"{body} centre", ec, 1e-12)
    check(f"{body} rotation", eR, 1e-12)
    check(f"{body} velocity and omega", ev, 1e-12)

    # Written surface at every write time against x = c + R (X0 - c0)
    eS = 0.0
    for row in d:
        tname = None
        for n in os.listdir(os.path.join(serial, "postProcessing",
                                         "rigidBodies")):
            try:
                if abs(float(n) - row[0]) < 1e-9*max(1, row[0]):
                    tname = n
            except ValueError:
                pass
        X = read_stl(os.path.join(serial, "postProcessing", "rigidBodies",
                                  tname, body + ".stl"))
        c, R, _, _ = path(row[0])
        if c is None:
            exact = X0
        else:
            exact = c + (X0 - c0) @ R.T
        if X.shape != exact.shape:
            fails.append(body + " STL size")
            continue
        eS = max(eS, np.abs(X - exact).max())
    check(f"{body} written STL", eS, 1e-7)

    if body == "hopperWall":
        check("hopperWall fixed (bit for bit pose)",
              float(np.abs(d[:, 1:13] - np.r_[d[0, 1:4], np.eye(3).ravel()]).max()),
              0.0)

    # 2 processors: master writes the same poses
    dp = read_dat(parallel, body)
    check(f"{body} 2 processors against serial",
          float(np.abs(dp - d).max()) if dp.shape == d.shape else 1.0, 0.0)

if fails:
    print("FAILED:", ", ".join(fails))
    sys.exit(1)
print("all rigid body checks passed")
