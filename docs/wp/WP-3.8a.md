# WP 3.8a: Axisymmetric wedge without the free surface

State: done
Branch: `claude/project-thread-53vkqx`

## Goal

Run a cylindrical problem as an OpenFOAM wedge slice (Sia 2026-11-18: Abaqus CEL needs
a quarter model): wedge patches in the Lagrange step, a time step that does not
collapse at the axis.

## Read first

- `docs/decisions/0061-axisymmetric-wedge.md`
- `docs/decisions/0018-lagrange-step-cpp.md` (boundary kinds), `0030-gravity-walls-cpp.md` (hydrostatic mirror)
- `src/celEulerian/eulerianSolid/eulerianSolidLagrange.C` (`kindOf`, the boundary
  face loop, `extrapolatedWalls`)

## Build

- `boundaryKind::wedge`: rotated cell state on wedge faces, no volume flux.
- `extrapolatedWalls` for the wall stress on wedge meshes.
- `eulerianSolid::cellSize()`: in-plane length with the radius for the time step.

## Done when

- `tests/test_wedge.sh` passes: Lame thick cylinder beyond two wall cells within 3 %
  (5 mm) and 1 % (2.5 mm), wall cells first order; lake at rest with axis cells below
  1e-10 m/s; cylindrical Janssen tank settles, K at the wall 0.27 to 0.45, mean s_v
  within 5 % of circular Janssen with that K (as the Cartesian Janssen check);
  2 processors match serial.

## Out of scope

- Free surface on a wedge (WP 3.8b), the silo on a wedge (3.8c), STL walls (3.8d).

## Log

- 2026-11-18: done. Four fixes on the way (decision 0061): the wedge face carries the
  predicted cell state (the old state sheared a falling plug at the axis); it keeps the
  full rotated velocity (dropping its normal part grew a swirl at the axis); wall
  stress extrapolated (hoop error at the wall was O(1)); the time step includes the
  radius (axis prisms were unstable at maxCo 0.4). The done-when row said Janssen
  within 2 %; the test uses 5 % like the Cartesian check (measured 2.0 %).
