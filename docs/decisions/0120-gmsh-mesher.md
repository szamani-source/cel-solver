# 0120: Gmsh as a second mesher in celSetup

Date: 2026-11-23. Thread "Build case generator and mesh checker" (WP G.1).

## Context

Sia (2026-11-23): meshing must be open source, and the setup tool "should handle
meshing either using gmsh or blockMesh". The validated Couto mesh is a four-block
wall-aligned blockMesh (decisions 0105, 0107, 0115) and must stay the default.

## Decision

- `mesh: mesher: blockMesh | gmsh` (default blockMesh). Both meshers take the same block
  layout from `conical.layout()`: block corners, the shell mid-line and domain side as
  (r, y) paths, column and row counts with their grading segments. The Gmsh path writes
  `system/mesh.geo` (built-in kernel): one transfinite, recombined quad surface per
  block and row segment (curves split where the grading changes, so each side has one
  `Progression`), splines for the curved block sides, then a one-layer revolution
  (wedge, half angle as blockMesh) or a one-layer extrusion (plane cases).
- Gmsh writes MSH 2.2; OpenFOAM's `gmshToFoam` converts it, and celSetup sets the
  patch types (wedge front/back, empty for plane) that gmshToFoam writes as `patch`.
  Physical groups name the same patches as blockMesh (outlet, bottom, outer, top,
  front, back), so materialProperties, the fill and celCheck do not change.
- Gmsh runs through its Python module when installed (`pip install gmsh`), else the
  `gmsh` command (`apt install gmsh`). Allrun meshes through `celSetup --mesh .`.
- Revolving a curve on the axis collapses it: Gmsh then makes prisms on the axis
  (158 for Couto) and drops that side from the extrusion result. The physical-group
  indices account for that.

## Consequences

- Couto with Gmsh: 7,584 cells (7,426 hex + 158 prisms on the axis), the same
  patches and face counts as blockMesh, checkMesh OK, non-orthogonality max 37.1 deg
  (blockMesh 34.7), celCheck time step 5.851e-6 s (blockMesh 5.850e-6 s), fill
  1503.4 kg, celFoam runs and keeps the mass. Interior nodes differ a little: Gmsh
  uses transfinite interpolation of the four sides, blockMesh its own; visible near
  the cylinder-hopper corner (`/mnt/project-files/setup/G1-blockmesh-vs-gmsh.png`).
- Gmsh is not in the CI image: `tests/test_celsetup_gmsh.sh` checks the written
  .geo there and skips the meshing, which runs wherever Gmsh is installed.
- The S.6 time-step options (`outletColumns: even`, `axisColumn`) are blockMesh only
  for now.
- The blockMesh mesh stays the validated one; a Gmsh discharge run is not part of
  G.1 (no long runs: accuracy first).
