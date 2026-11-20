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
