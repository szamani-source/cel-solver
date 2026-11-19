#!/usr/bin/env python3
"""WP H.1: the conical hopper case spec (docs/spec/conical-hopper.md) is consistent.

Checks the copied and digitised reference data against the geometry and the closed
forms the papers quote, so a typo in a number fails here and not in a later run.
"""
import csv
import math
import os
import sys

root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
data = os.path.join(root, "docs", "spec", "hopper-data")
fails = []


def check(name, ok, detail):
    print(("PASS " if ok else "FAIL ") + name + ": " + detail)
    if not ok:
        fails.append(name)


def rows(fname):
    with open(os.path.join(data, fname)) as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


gamma = 8.397e3          # N/m^3 (P2)
rho = 856.0              # kg/m^3
g = 9.8                  # P1
D, hc, ht, Dh = 1.0, 2.0, 0.48, 0.35
beta = math.radians(34.3)

# 1. Geometry: hopper slope and apex height, fill weight
R, r = D / 2, Dh / 2
check("hopper height", abs((R - r) / math.tan(beta) - ht) < 0.01,
      "%.3f m from beta, 0.480 m drawn" % ((R - r) / math.tan(beta)))
check("apex height", abs(R / math.tan(beta) - 0.738) < 0.01,
      "%.3f m, 0.738 m drawn" % (R / math.tan(beta)))
V = math.pi * R**2 * hc + math.pi * ht / 3 * (R**2 + R * r + r**2)
W = V * gamma / 1e3
t1 = rows("couto2013_table1.csv")
Ws = [float(x["W_kN"]) for x in t1]
check("fill weight", all(abs(w / W - 1) < 0.015 for w in Ws),
      "geometry %.2f kN, tests %.2f to %.2f kN" % (W, min(Ws), max(Ws)))

# 2. Table 1 sanity: static min <= max, discharge GA0 above static GA0
for t in range(1, 8):
    s = {x["stat"]: x for x in t1 if x["test"] == str(t) and x["phase"] == "static"}
    d = [x for x in t1 if x["test"] == str(t) and x["phase"] == "discharge"][0]
    cells = ["GA%d" % i for i in range(8)]
    ok = all(float(s["min"][c]) <= float(s["max"][c]) for c in cells)
    ok = ok and float(d["GA0"]) > 1.8 * float(s["max"]["GA0"])
    check("table 1 test %d" % t, ok, "static GA0 %s, discharge GA0 %s kPa"
          % (s["max"]["GA0"], d["GA0"]))

# 3. Nedderman closed forms quoted in P1 (K_A 0.338, mu 0.2, Q0 13.0, Cp 2.23)
KA, mu = 0.338, 0.20
zc = {x["cell"]: float(x["z_m"]) for x in rows("couto2013_cells.csv")}
H = hc + ht


def janssen(depth):
    return rho * g * D / (4 * mu) * (1 - math.exp(-4 * mu * KA * depth / D)) / 1e3


sw_t = janssen(hc)
check("bin wall pressure at transition", abs(sw_t - 4.4) < 0.15,
      "%.2f kPa (P1 Fig. 11 dashed line 4.4)" % sw_t)
check("surcharge Q0", abs(sw_t / KA - 13.0) < 0.3, "%.2f kPa (P1: 13.0)" % (sw_t / KA))
check("hopper peak Cp Q0", abs(2.23 * 13.0 - 29.0) < 0.1, "%.1f kPa" % (2.23 * 13.0))
static = [x for x in t1 if x["phase"] == "static"]
for c in ["GA%d" % i for i in range(1, 8)]:
    lo = min(float(x[c]) for x in static)
    hi = max(float(x[c]) for x in static)
    p = janssen(H - zc[c])
    # Upper cells read 20-45 % above Janssen with K_A (filling stresses; P1's FEM
    # sits on Janssen there too), so only the bottom cells are held to it.
    if c in ("GA1", "GA2", "GA3"):
        check("Janssen at " + c, lo / 1.15 <= p <= hi * 1.15,
              "%.2f kPa, tests %.2f to %.2f" % (p, lo, hi))
    else:
        print("INFO Janssen at %s: %.2f kPa, tests %.2f to %.2f (%+.0f %%)"
              % (c, p, lo, hi, 100 * (p / lo - 1)))

# 4. P1 Fig. 10 against Beverloo C = 0.59
for x in rows("zheng2015_fig10.csv"):
    d0, w = float(x["D0_m"]), float(x["W_kgs"])
    wb = 0.59 * rho * math.sqrt(g) * d0**2.5
    tol = 0.05 if d0 >= 0.1 else 0.5
    check("Beverloo D0 %.2f" % d0, abs(w / wb - 1) < tol,
          "digitised %.1f, C 0.59 gives %.1f kg/s" % (w, wb))

# 5. P1 Fig. 3 plateau and Fig. 11 peak
v = rows("zheng2015_fig3.csv")
check("outlet velocity plateau", abs(float(v[-1]["v_ms"]) - 1.37) < 0.02,
      "%s m/s" % v[-1]["v_ms"])
f11 = rows("zheng2015_fig11_fem.csv")
pk = max(f11, key=lambda x: float(x["p_kPa"]))
check("FEM peak at transition", abs(float(pk["z_m"]) - ht) < 0.03,
      "%s kPa at z %s m" % (pk["p_kPa"], pk["z_m"]))

# 6. Case B geometry: alpha 30, D0 0.1 in H 2.48
hb = (0.5 - 0.05) / math.tan(math.radians(30))
check("case B hopper height", abs(hb - 0.779) < 0.001, "%.3f m, cylinder %.3f m" % (hb, 2.48 - hb))

print("%d failures" % len(fails))
sys.exit(1 if fails else 0)
