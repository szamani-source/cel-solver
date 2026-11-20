# 0098: The Taylor-bar lip is first-order smearing; check the extrapolated widening

Date: 2026-11-19 (accuracy package A.7, thread "Fix sand footing bearing capacity")
Status: accepted (test check added; no solver change).

## Context

Decision 0063: the Eulerian Taylor bar on the wedge reads the widening 16 % short at
0.2 mm (5.02 against the Abaqus CAX4R 6.006 mm), with the thin lip at the impact face
smeared over a cell or two.

## Findings (tests/test_taylor_wedge.sh, 2026-11-19)

| Cells | Widening, bottom row (mm) | Bottom row extrapolated to the floor | Shortening (mm) |
| --- | --- | --- | --- |
| 0.4 mm | 3.66 | 3.73 | 12.47 |
| 0.2 mm | 5.01 | 5.10 | 13.02 |
| h -> 0 (Richardson, first order) | 6.37 (+6.0 %) | | 13.57 (+3.5 %) |

The radius falls steadily with height above the floor (5.01, 4.83, 4.71, 4.56 mm in the
first four rows at 0.2 mm). Measuring at the floor instead of the row centre adds only
2 %, so measuring at the floor is not the cause. The error halves with the cell size, as
a volume fraction smearing a lip that is thin compared with a cell would. The converged
Eulerian widening is within about 6 % of the Lagrangian reference, whose lip is itself
mesh-sensitive: the same benchmark's 3D C3D8R run gives 5.528 mm (validation row 11).

## Decision

1. No solver change: a sharper lip needs interface reconstruction (PLIC) in the
   remap, which is not planned.
2. `tests/taylorwedge/taylor_wedge_checks.py` adds a check: the widening extrapolated
   to zero cell size from the 0.4 and 0.2 mm runs within 10 % of 6.006 mm. The 20 %
   check at 0.2 mm stays.

## Consequences

- Thin flanges and lips need several cells across them, or the Lagrangian body.
- Known issue: the Taylor-bar lip is closed as a resolution effect.
