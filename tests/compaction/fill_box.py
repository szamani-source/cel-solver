"""WP K.3 study (decision 0113): does the cap remove the doubled at-rest
pressures of a box filled from the top (decision 0108)? Not a pass/fail test.

The WP 3.7a inlet box (tests/inlet/inlet_checks.py), 0.6 m tall, fed at
2 m/s through the 40 mm lid inlet for 1 s (T0 rules plus fullBlend), damping
from 1.4 s, end 2.4 s; against the same mass placed at once with gravity
ramped over 0.5 s. With and without the cap (lambda* = 0.02, LAM env).
Prints the lateral ratio K = s_xx/s_yy in the middle of the deposit and the
right-wall normal stress below mid-height.

    python3 fill_box.py WORKDIR gravity|top cap|nocap
"""
import os, sys, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "inlet"))
import inlet_checks as ic
vc = ic.vc
work, mode, cap = sys.argv[1], sys.argv[2], sys.argv[3] == "cap"
LAM = float(os.environ.get("LAM", "0.02"))
# step 1 of the top-fill fix (Sia 2026-11-21): slower landing. S_IN sets the
# inlet speed (feed time scales so the mass is the same), BOX_H the box height
# (drop height)
S_IN = float(os.environ.get("S_IN", "2.0"))
T_FEED = 2.0 / S_IN
ic.H = float(os.environ.get("BOX_H", "0.6")); ic.S_IN = S_IN; ic.T_START = 0.05; ic.T_STOP = 0.05 + T_FEED
T_DAMP, T_END = T_FEED + 0.35, T_FEED + 1.35
name = f"{mode}_{'cap' if cap else 'nocap'}_s{S_IN:g}_h{ic.H:g}" + (f"_relax{os.environ['RELAX']}" if os.environ.get("RELAX") else "")
p = ic.setup(work, name)
mp = f"{p}/constant/materialProperties"; t = open(mp).read()
extra = "fullBlend yes;\n" + (f"capLambda {LAM!r};\ncapP0 100;\n" if cap else "")
# step 2 (decision 0114): RELAX = rest relaxation time [s] (speed 0.01 m/s)
RELAX = float(os.environ.get("RELAX", "0"))
if RELAX > 0:
    extra += f"restRelax {{ time {RELAX!r}; speed 0.01; }}\n"
t = t.replace("normalCap yes;\n", "normalCap yes;\n" + extra)
Hf = 0.04 * 2.0 * 1.0 / 0.2      # same mass as the fed box
if mode == "gravity":
    t = t.replace(f"endTime {ic.T_STOP!r};", "endTime 0.0;").replace("startTime 0.05;", "startTime 0.0;")
open(mp, "w").write(t)
dct = f"{p}/system/voidCheckDict"; d = open(dct).read()
d = d.replace("tEnd 1.2;", f"tEnd {T_END};")
import re
d = re.sub(r"tEnd [^;]*;", f"tEnd {T_END};", d)
if mode == "gravity":
    d = re.sub(r"fill \([^)]*\);", f"fill (0 0 {ic.W} {Hf});", d)
open(dct, "w").write(d)
open(f"{p}/system/controlDict", "a").write(f"damping 100;\ndampingStart {T_DAMP if mode=='top' else 0.0};\ngravityRamp {0.0 if mode=='top' else 0.5};\n")
vc.run(p)
c = np.loadtxt(f"{p}/cells.dat")
x, y, al, rho = c[:,0], c[:,1], c[:,2], c[:,3]
sxx, syy = c[:,6], c[:,7]
h = np.loadtxt(f"{p}/history.dat")
full = al > 0.99
top = max(y[full]) if full.any() else 0
mid = full & (y > 0.2*top) & (y < 0.7*top) & (x > 0.06) & (x < 0.14)
K = (sxx[mid]/syy[mid]).mean()
wall = full & (x > ic.W - 0.01) & (y < 0.5*top)
print(f"{name}: mass {h[-1,2]:.3f} kg, top {top:.3f} m, K mid {K:.2f}, wall sxx at y<{0.5*top:.2f}: "
      + " ".join(f"{v:.0f}" for v in sxx[wall][::2]) + f"; mean rho {np.sum(rho*al)/np.sum(al):.0f}, Umax {np.abs(c[full,4:6]).max():.1e}", flush=True)
