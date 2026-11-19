# 0103: Angled (non-orthogonal) hex cells: time step and plan

Date: 2026-11-19 (WP G.1)

## Context

Sia (2026-11-19) wants the hopper meshed with wall-aligned blockMesh blocks
(trapezoid blocks, cells parallel to the hopper wall) instead of a Cartesian mesh
cut by an STL wall, and said it is very important that the solver handles angled
hex cells of varying shape in general. The audit
(`docs/spec/angled-cells-audit.md`) found that the core scheme already works on
angled cells, but the time step and every free-surface, void, wall and STL rule
assume faces along x, y or z.

## Decision

1. The solver detects angled meshes at start: a face (empty and wedge sides
   excepted) whose normal is off the x, y and z axes, on any processor.
2. On an angled mesh the stable length of a cell is V over its largest face area
   (empty and wedge sides excepted): the perpendicular height between the
   closest pair of opposite faces, exact for a parallelepiped. cbrt(V) (the
   Cartesian length) overestimates it on sheared or thin cells: by sqrt(2) on a
   45 degree parallelogram of square footprint. On a wedge, this in-plane length
   is combined with the radius for the hoop strain as before (decision 0061).
3. Axis-aligned meshes keep the old length, so every existing run is unchanged.
4. The work is split into packages G.1 to G.6 (`docs/work-packages.md`); every
   one keeps Cartesian meshes bit-for-bit and replaces x/y/z directions with each
   cell's face-pair (mesh-line) directions, which reduce to x/y/z on a Cartesian
   mesh.

## Consequences

- `tests/test_angled_cells.sh` (CI): P-wave on a 45 degree sheared mesh within
  5 % of the Cartesian error, second order; elastic column on leaning cells of
  varying shape with the solver's own step settles to the closed-form stress
  within about 1 %.
- freeSurface still stops with a clear error on angled meshes until G.2.
