"""WP 4.0: assert the recommended coupling (cut-cell moving-wall Riemann) on the
prototype experiments of prototype/contact2d.py, and that the penalty coupling it
was compared with is still the noisier one (decision 0042)."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "prototype"))
import contact2d as c  # noqa: E402

fails = []


def check(name, val, lo, hi):
    ok = lo <= val <= hi
    print(f"{'ok  ' if ok else 'FAIL'} {name}: {val:.4g} (band {lo:g} to {hi:g})")
    if not ok:
        fails.append(name)


for off in (0.0, 0.5):
    r = c.piston("riemann", offset=off)
    check(f"riemann piston offset {off}: stress / z v", r["mean"], 0.98, 1.02)
    check(f"riemann piston offset {off}: noise", r["noise"], 0.0, 0.01)
    check(f"riemann piston offset {off}: leak", r["leak"], 0.0, 1e-3)
    check(f"riemann piston offset {off}: mass error", r["mass"], 0.0, 1e-12)
    rn = r["noise"]
for off in (0.0, 0.5):
    r = c.block("riemann", offset=off)
    check(f"riemann block offset {off}: leak", r["leak"], 0.0, 1e-3)
    check(f"riemann block offset {off}: mass error", r["mass"], 0.0, 1e-12)
r = c.incline("riemann")
check("riemann incline 30 deg, delta 20: acceleration / exact", r["acc"], 0.97, 1.03)
check("riemann incline: leak", r["leak"], 0.0, 1e-4)

p = c.piston("penalty", offset=0.0, pfac=1.0)
print(f"penalty k = M/h piston: stress / z v {p['mean']:.4f}, noise {p['noise']:.4f},"
      f" leak {p['leak']:.2e}")
check("penalty noise over riemann noise", p["noise"] / max(rn, 1e-12), 10.0, 1e9)

if fails:
    print("FAILED:", ", ".join(fails))
    sys.exit(1)
print("all contact scheme checks passed")
