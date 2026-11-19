# 0063: Eulerian Taylor bar on the wedge, widening tolerance 20 %

Date: 2026-11-18 (WP 2.6)

## Context

The Taylor bar (copper rod 32.4 mm long, radius 3.2 mm, 227 m/s onto a rigid
frictionless floor, von Mises 314 MPa perfectly plastic, 80 microseconds) is
axisymmetric, so it runs on the wedge (decision 0061) with void around the rod. The
reference is the axisymmetric CAX4R run of the Abaqus Benchmarks ("Impact of a copper
rod"): shortening 13.11 mm, widening of the impact face edge 6.006 mm. The Lagrangian
quarter model (WP 5.2, decision 0071) is checked at 5 % and 10 %.

| cells | shortening (mm) | widening (mm) | run |
|---|---|---|---|
| 0.4 mm | 12.47 (-4.9 %) | 3.67 (-39 %) | 45 s |
| 0.2 mm | 13.02 (-0.7 %) | 5.02 (-16 %) | 6.5 min |

Final length from the material in the axis column; face radius from the bottom cell
row, area weighted. Mass kept to 1e-14.

## Decision

`tests/test_taylor_wedge.sh` runs both sizes and checks shortening within 5 % and
widening within 20 % at 0.2 mm, both errors falling with the cell size, and the bar
at rest (kinetic energy below 5 %).

The widening is short because the face edge spreads into a thin lip that a volume
fraction smears over a cell or two, and the bottom row sits half a cell above the
face. It roughly halves its error when the cells halve (first order), so 0.1 mm should
come near 10 %; at about 50 min a run that is left out of CI.

## Consequences

- Eulerian results for thin flanges and lips need several cells across them; the
  Lagrangian body (WP 5.x) is the better tool for those.
- No hardening, rate or temperature terms (Johnson-Cook) yet; the benchmark is
  perfectly plastic.
