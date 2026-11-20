# 0107: STL walls with a free surface on angled cells (WP G.5)

## Context

Sia (2026-11-20) asked for the Couto silo with the wall as an STL shell
inside a wall-aligned mesh, with a few cells of buffer outside it, so that a
leak through the wall shows up as material in the buffer. The solver
stopped with an error for STL walls plus a free surface on angled cells.

## Decision

- The error is removed. The one per-axis rule left in the cut-cell set-up,
  the open length used by the remap's slab rule (open volume over the larger
  open face area across an axis), now gives each face's area to the mesh
  line it lies across (`inLines`, as decision 0104 did for the free
  surface). The wall face fluxes already used `inLines` (G.2).
- Axis-aligned meshes take no new path: their results are unchanged.

## Test

`tests/test_angled_stl.sh`: lakes at rest on leaning cells against the
Cartesian mesh, in a tank of vertical STL walls and against an STL wall at
30 degrees: mass to 1e-14, nothing past the walls, flat surface, s_yy within
0.6 % of hydrostatic, wall-normal stress within 0.9 % of the Cartesian run.

## Consequences

- `tests/hopper/buffered_case.py` builds the shell (20 mm, fillets of 50 mm
  at the transition and 15 mm at the outlet edge, a 20 mm spout) inside a
  wall-aligned mesh (`aligned_buffered_case`) or a graded Cartesian one
  (`buffered_case`), with a leak report; `tests/manual/test_hopper_pressure.sh
  mc aligned_buffered` runs it.
- Only lakes at rest are tested here; the silo run is the first flowing case.
