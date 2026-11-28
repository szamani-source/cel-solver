# 0140: Tutorials in three category folders; solid-only set

Date: 2026-11-27 (WP TU.1)

## Context

Sia (2026-11-26) asked for tutorial cases built from the validated cases, as templates
for users, and a professional user manual organised in three categories: solid only,
solid-fluid, and solid-fluid-temperature. The tutorials had grown one folder per work
package at the top of `tutorials/`, with READMEs of different shapes and some results
from older versions of the solver.

## Decision

- `tutorials/` has one folder per category: `solid/`, `solidFluid/`,
  `solidFluidThermal/`. `tutorials/README.md` is the index. The solid tutorials moved
  into `tutorials/solid/` under their old names, so earlier notes still find them by
  name; tests, the installer, `celfoam tutorials` and the manual follow the new paths.
- The solid set is `smokeBox`, `siloDischarge`, `siloCohesion`, `siloFed`, `stlSilo`,
  `stlHopper`, `coutoHopper`, `setupPlaneTank`, `stripFooting` (new: the WP 4.4
  Prandtl footing, the coupled Eulerian-Lagrangian demonstration) and
  `rigidBodyMotion`. The old `stlTank` (a 0.3 s settling demo with no README) is
  removed: `stlSilo` covers it.
- Every tutorial README has the sections Setup, Numerics, How to run and Expected
  results, and shows at least one result figure, so the manual can take them as they
  are. `tests/test_solid_tutorials.sh` checks this and runs a short stretch of each case.
- Figures come from full runs of the tutorials as shipped, on the current solver. The
  Couto hopper is the exception: its figures are from the validated 4 mm runs of
  2026-11-23 (PR #78, decision 0117), which take 8 h on 4 cores; nothing in the solid
  solver has changed since that affects it.
- The compaction cap gets no tutorial: it is checked at one point and in a box, but it
  did not close the at-rest gap in the Couto silo (Sia closed the item 2026-11-24).
- `hotSiloDischarge` stays at the top level until the solid-fluid-temperature tutorials
  are gathered, because the temperature thread (T.7) is still changing it.

## Consequences

- Installed copies put the tutorials in `<dest>/solid/<name>`; the first command in
  README.md, INSTALL.md and the manual's quick start changed accordingly.
- The gas and temperature tutorials go into `solidFluid/` and `solidFluidThermal/` with
  the same README sections; a test like `tests/test_solid_tutorials.sh` can cover them.
