# 0062: Wedge silo against a quarter model needs three cells across the outlet

Date: 2026-11-18 (WP 3.8c)

## Context

WP 3.8c compares a cylindrical silo (R = 0.1 m, fill 0.2 m, T0 Mohr-Coulomb with the
T0 rule set, 10 mm cells) run on the axisymmetric wedge (decision 0061) with the same
tank as a 3D quarter model (symmetry planes, STL cylinder wall, square quarter
outlet). The done-when asks for the two within 10 % and both in the 3D Beverloo band.

With the first outlet, D = 0.04 m (2 cells across the radius, a 2 x 2 cell quarter
outlet), the wedge read C = 0.688 (0.667 at 5 mm) and the quarter C = 0.495, a 28 %
gap. Both full tanks matched (11.30 kg).

With D = 0.06 m (3 cells; quarter outlet 0.03 x 0.03 m, compared by the circle of
equal area, D_eq = 0.0677 m) the wedge reads C = 0.678 and the quarter 0.635: 6.8 %.

## Decision

The test uses the 0.06 m outlet. The gap at 2 cells is put down to the quarter
outlet being under-resolved: a 2 x 2 square loses a larger share of its area to the
slow edge cells than the wedge's 2-cell radius, which keeps its full annular area.
The wedge result barely moves with the outlet size or the cell size, so it is the
better-resolved of the two.

## Consequences

- Quarter (or full 3D) silo models need at least three cells across the outlet
  half-width; the wedge needs only two.
- The wedge ran 33 s against 164 s for the quarter on one processor (5x) at 10 mm.
- Not checked: a 5 mm quarter model (about 16x the cost), and non-central outlets,
  which a wedge cannot represent.
