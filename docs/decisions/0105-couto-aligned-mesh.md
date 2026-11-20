# 0105: Couto silo on the graded wall-aligned mesh

Date: 2026-11-20 (WP H.7)

## Context

Sia asked for a wall-aligned hopper mesh (decision 0103) and then for cells finer
at the outlet and the walls, coarse at the top, to match the papers. On the uniform
aligned mesh (25 mm) the discharge crept at 21 kg/s (8 kg/s at 50 mm). Controls ruled
out wall friction, flow impedance, the outlet boundary (a free-fall zone below the
outlet did not help: the hopper itself held back) and dilation (psi 1, 5, 15 gave the
same rate).

## Findings (25 mm base cells, wedge, both models with the T0 rules)

| Mesh and rules | Rate | Largest speed | GA0 discharge |
| --- | --- | --- | --- |
| Uniform aligned | 21 kg/s | 0.4 m/s | 21.3 kPa |
| Outlet rows 3x finer | 63-86 kg/s | 20-66 m/s, stress to 4.7 MPa | |
| Outlet and walls 3x finer | 74-99 kg/s | 23-52 m/s | 18.4 (3-34) |
| Outlet rows 3x finer, fullBlend | 85 kg/s | 1.4 m/s | 19.7 |
| Same, walls also 3x finer | 84 kg/s | 1.4 m/s | 19.7 |
| Same, DP compression fit | 45 kg/s | 0.8 m/s | 22.4 |
| Same, DP plane-strain fit | 86 kg/s | 1.4 m/s | 19.4 |

The outlet rows set the rate. The spikes on the graded mesh start as compression
of full cells beside loosened ones (density +3 %); `fullBlend yes` (decision 0049)
removes them. Beverloo as P1 uses it (C 0.59, k 0) gives 115 kg/s.

## Decision

- `tests/manual/test_hopper_pressure.sh [mc|dp] [aligned|cut]` runs the aligned mesh
  by default: `hopper_case.aligned_case`, outlet rows 3x finer, `fullBlend yes`
  added to the decision 0102 rules. `cut` keeps the H.3 set-up.
- Drucker-Prager is run with the plane-strain fit for this case.

## Results to 1.5 s (`CEL_NPROCS=4`)

- MC: 84.7 kg/s, largest speed 1.4 m/s. At rest GA1-GA5, GA7 pass; GA0 11.7 kPa
  (+21 %), GA6 -14 %. Discharge GA0 19.2 and GA1 6.8 kPa pass; GA2 4.0 (-32 %).
- DP (plane strain): 86.9 kg/s, 1.4 m/s. At rest GA2-GA7 pass, GA1 +23 %, GA0 +17 %.
  Discharge GA0 19.0, GA1 6.7 pass; GA2 4.8 (-18 %).

## Open

GA0 at rest (same on every mesh: the elastic fill's stress ratio), GA2 in discharge,
and the rate gap to Beverloo / Sia's Abaqus runs.
