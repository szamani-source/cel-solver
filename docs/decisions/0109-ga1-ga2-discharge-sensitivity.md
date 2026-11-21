# 0109: GA1/GA2 low in discharge is not a model knob

## Context

In discharge the Couto cylinder cells just above the transition read low on
every mesh: GA1 5.3-6.8 kPa against 7.1-8.2, GA2 3.9-4.6 against 5.8-6.5
(fine aligned, rounded-corner STL and 25 mm aligned meshes, MC and DP, both
starting states). Hopper wall cells move at 0.2-0.3 of the axis speed (radial
flow in this cone gives about 0.7).

## Sensitivity (MC, 25 mm aligned mesh, outlet rows 3x finer, T0 rules,
discharge mean 0.25-1.0 s after opening)

| change | rate kg/s | GA0 | GA1 | GA2 |
|---|---|---|---|---|
| reference (wall mu 0.2) | 84.9 | 19.2 | 6.69 | 3.89 |
| wall mu 0.1 | 88.0 | 26.0 | 6.44 | 4.21 |
| wall mu 0.05 | 93.9 | 28.4 | 6.34 | 4.38 |
| without normalCap | 57.4 | 20.1 | 6.81 | 4.05 |
| without flowImpedance | 64.6 | 19.7 | 6.53 | 3.98 |
| without shearCap | 84.3 | 19.1 | 6.66 | 3.94 |
| dilation psi 5 deg (from 1) | 60.0 | 20.2 | 6.78 | 4.55 |
| starting gap (rounded STL mesh) | 122.9 | 20.1 | 5.77 | 4.32 |
| Couto | (110-131 Beverloo) | 18.7-24.8 | 7.1-8.2 | 5.8-6.5 |

Wall friction moves GA0 and the rate, not GA1/GA2. No T0 rule, dilation or
starting state moves GA1/GA2 by more than 0.5 kPa.

## Decision

Stop tuning for GA1/GA2. The cylinder above the transition stays an elastic
plug in every run; the switch pressure needs the converging flow zone to
reach above the transition, which a shear band width set by the mesh
cannot give. Revisit after the mesh-dependence package (law with a length
scale), which Sia scheduled before the Linux install.

## Consequences

GA1/GA2 stay a known gap in STATUS. Runs: /tmp scratch (wallmu.py,
rules.py variants of hopper_pressure_checks), about 18 min each on 4 cores.

## Revisit with the nonlocal law (2026-11-21)

Best setup (rounded-corner STL mesh, outlet cells 3-5 mm, 25 mm in the
upper cylinder, 2-cell wall-gap start, `gap_case.py`), MC, the T0 rules
plus fullBlend, with the nonlocal fluidity of decision 0111 (A 0.48,
d 4 mm). Discharge mean 0.25-1.0 s after opening, kPa:

| run | rate kg/s | GA0 | GA1 | GA2 | GA3 |
|---|---|---|---|---|---|
| local mu(I) (above) | 122.9 | 20.1 | 5.77 | 4.32 | |
| nonlocal, spreadLayers yes | 119.4 | 24.1 | 5.63 | 5.03 | 4.51 |
| nonlocal, spreadLayers no | 119.3 | 24.3 | 5.31 | 5.18 | 5.27 |
| Couto | (110-131 Beverloo) | 18.7-24.8 | 7.1-8.2 | 5.8-6.5 | 4.4-5.6 |

GA2 moves up about 0.8 kPa (now 11-14 % low), GA1 does not move (25-30 %
low), the rate stays inside Beverloo, no leak. About 2 h each on 4 cores.
The cooperative length is a few grain sizes (about 3 d = 12 mm), below the
25 mm cylinder cells, so the law cannot widen the converging zone there;
resolving it would need cells near d in the cylinder (about 25 times more
cells). GA1 stays a known gap; it is not a reason to hold the Linux
install. Abaqus (Zheng and Yu) is also low at GA2.
