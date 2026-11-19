# WP 3.8b: Free surface on the wedge

State: done
Branch: `claude/project-thread-53vkqx`

## Goal

`freeSurface yes` (void, outlets, inlets) on an axisymmetric wedge mesh, so a
cylindrical silo can be filled, settle and later discharge (3.8c).

## Read first

- `docs/decisions/0061-axisymmetric-wedge.md` (with its 3.8b addendum)
- `docs/decisions/0031-free-surface-cpp.md` (slab rules, void closure)

## Build

- Wedge sides pass the Cartesian check in `setCartesianGeometry` (hoop length).
- `hoopRate_`: the u/r part of the radial flux goes to the hoop axis in the void
  closure and in the remap's slab lengths.

## Done when

- `tests/test_wedge_void.sh` passes: lake at rest (1e-8 m/s, mass 1e-12); falling
  plug at g within 0.5 %, lands, mass 1e-12; T0 void-top tank within 5 % of circular
  Janssen, holds when undamped; 2 processors match serial.

## Out of scope

- Outlet on the axis and the silo discharge (3.8c); STL walls on a wedge (3.8d).

## Log

- 2026-11-18: done. A free-falling plug sheared at the axis until the wedge faces
  carried the predicted cell state (fixed in 3.8a's code). The void-top tank needs the
  T0 rule set; with plain MC it never settles, on the Cartesian mesh too. An elastic
  plug's landing is chaotic at the axis (waves focus there), so the parallel check
  uses the tank. Done-when row said "as `tests/test_void.sh`"; the translating and
  separating blocks were left out (Cartesian-only motions, nothing wedge-specific).
