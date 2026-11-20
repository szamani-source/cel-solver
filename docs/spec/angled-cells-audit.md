# Angled (non-orthogonal) hex cells: audit and plan

Sia (2026-11-19): hopper meshes will be hex cells of varying shape, aligned with
the walls (blockMesh trapezoid blocks). The solver must handle angled cells.

## What already works (probes, 2026-11-19)

- Core FV scheme (stress divergence, velocity update, remap without void): a
  P-wave on a 45 deg sheared periodic mesh converges at 2nd order with the same
  error as the Cartesian mesh (L1 1.65e-2, 3.84e-3, 8.17e-4 at n 25, 50, 100).
- Couto silo at rest on a wall-aligned two-block wedge mesh (lid, no free
  surface): settles to 3e-4 m/s, cylinder pressures within 1 to 2 % of Janssen
  (K 0.43), smooth hopper wall pressure 10.1 to 11.9 kPa.

## Cartesian assumptions found

| Where | Assumption | Effect on angled cells | Package |
|---|---|---|---|
| `cellSize()` (eulerianSolid.C) | length = cbrt(V) off the wedge | time step too large on thin or sheared cells | G.1 |
| `setCartesianGeometry` (eulerianSolid.C:728) | every face normal is a coordinate axis | fatal error with freeSurface | G.2 |
| `faceAxis()`, `slabSides()` (eulerianSolidLagrange.C:940) | face pairs per x/y/z axis | free surface, void closure, slab sides | G.2 |
| void closure dVa (eulerianSolidLagrange.C ~2528-2669) | swept volume along an axis | leakage at angled surfaces | G.2 |
| remap `axisOf` (eulerianSolidRemap.C:914-1163) | swept slices and H per axis | wrong partial-cell fluxes | G.2 |
| spreadLayers gAxis (eulerianSolidLagrange.C:1678) | gravity along a mesh axis | layers on angled columns | G.3 |
| extrapolatedWalls, equilibriumJump | faces normal to or along gravity | wall and surface stresses | G.3 |
| wall friction, outlet, inlet patches | patch normal is an axis | friction on angled hopper wall | G.3 |
| wedge `axisOf` (eulerianSolid.C:369) | trapezoid block fans to the axis | axis cells and hoop length | G.4 |
| immersed (STL) open lengths | per axis on a Cartesian mesh | STL on a skewed mesh | G.5 |
| MPI, restart | untested on skewed meshes | | G.6 |

## Rule

Every change keeps Cartesian meshes bit-for-bit (all existing tests unchanged).
Directions come from each cell's opposite face pairs (mesh lines), which reduce
to the x/y/z axes on a Cartesian mesh.

## Progress

- G.1 done (decision 0103): stable step from V over the largest face.
- G.2 done (decision 0104): mesh-line frames per cell; `cellSize`,
  `setCartesianGeometry`, `slabSides`, `touches`, void closure, remap and
  rule 9 use them.
