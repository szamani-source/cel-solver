# 0041: STL walls and bodies as immersed rigid surfaces on a Cartesian mesh

Date: 2026-11-17 (WP 4.1)

## Context

Sia (2026-11-17) asked for wall interactions defined by STL geometry (for example silo
walls and hopper inserts that do not follow the mesh), and whether the mesh must be hex
or could be polyhedral or hex-dominant (snappyHexMesh).

What the code assumes about the mesh today (checked in WP 4.1):

- `freeSurface yes` (void, outlet, shear cap, rule 9, `normalCap`; everything the silo
  needs) stops at start-up unless every face is normal to an axis
  (`eulerianSolid::setCartesianGeometry`, decision 0031). The void is a slab per axis,
  not isoAdvector. Graded blockMesh hexes pass; snapped cells and polyhedra do not.
- Without void, the Lagrange step (Gauss gradients, face Riemann solver along the face
  normal) and the remap (least-squares gradients, face fluxes) are written face by face
  and run on any OpenFOAM mesh, but were only tested on uniform hexes. Two places lose
  accuracy on other meshes: face states are extrapolated with half the centre-to-centre
  vector (assumes the face midway between the centres), and the stable time step uses
  cbrt(V), too large for flat cells.

Two ways to bring STL walls in: (a) mesh them conformally with snappyHexMesh and keep
walls as mesh boundaries, which needs general-mesh support first (PLIC void on general
cells, void closure along face normals, the two fixes above: about 3 to 4 packages);
(b) keep the Cartesian Eulerian mesh and let the STL surfaces cut through it, as
Abaqus/Explicit CEL and LS-DYNA ALE do (structured Eulerian mesh, Lagrangian walls).

## Decision

(b). STL walls and moving bodies are immersed rigid surfaces (library `celLagrangian`,
class `rigidBody`), coupled to the material by the scheme WP 4.0 recommends. The
Eulerian mesh stays Cartesian hex. Mesh-boundary walls (decision 0030) remain for walls
that follow the box.

WP 4.1 builds only input and kinematics:

- `constant/rigidBodies`, `bodies { <name> { surface <file>; centre; velocity; omega; } }`,
  surfaces from `constant/triSurface` (any format OpenFOAM's `triSurface` reads).
- Motion is prescribed by `Function1<vector>` velocity of a reference point and angular
  velocity. Position is the exact integral of velocity; rotation is the Rodrigues
  rotation by the integral of omega over the step, exact while the axis is fixed
  (second order otherwise). No forces yet.
- Every processor holds the whole surface (walls and inserts are small); only the master
  writes `postProcessing/rigidBodies/<name>.dat` and the moved STL at write times.
- A run that starts later than time zero is moved from the STL pose in one step (exact
  for the same cases). Restart of a body with a varying rotation axis is not supported.

## Consequences

- Polyhedral and hex-dominant meshes stay unsupported for the silo physics. If they are
  ever wanted, it is the package list in Context (a), recorded under Phase 3B.
- The coupling must handle a wall that cuts cells at any angle: WP 4.0 adds a fixed
  inclined frictional wall to its experiments, and WP 4.3a checks Janssen in a tank whose
  STL walls are tilted against the mesh.
- STL single precision (about 1e-8 m at 1 m) is the geometric accuracy of a body.
- Open surfaces (one plate as a wall) have no inside; WP 4.2 decides whether walls must be
  closed volumes or one-sided surfaces with a thickness.
