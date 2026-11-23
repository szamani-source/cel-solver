# 0120: User-provided meshes and boundary conditions in celSetup

Date: 2026-11-23. Thread "Build case generator and mesh checker" (WP G.3).

## Context

Sia (2026-11-23): "the GUI should have an option for loading a completely different
mesh (user provides) as well, where we can then choose the boundaries and assign
boundary conditions."

## Decision

- A template `userMesh` (`tools/celsetup/usermesh.py`). `mesh: file` is an OpenFOAM
  case or polyMesh folder, or a file OpenFOAM converts: .msh (Gmsh, told from Fluent by
  its `$MeshFormat` header), .unv, .neu, .ans; a .geo file is meshed with Gmsh first.
  `scale` runs transformPoints (meshes drawn in mm). The import is the case's mesh step
  (`celSetup --mesh`, Allrun), so the case folder holds its own copy.
- `boundaries:` gives every patch one condition: wall (friction), outlet (friction,
  openTime), inlet (speed, start/end, density, pressure, friction), symmetry
  (symmetryPlane), wedge, empty, free (plain traction-free patch). They map onto the
  existing materialProperties walls/outlets/inlets and patch types: no solver change.
  Every patch must have one (an unassigned patch is an error that names it); so must
  every named patch exist.
- The GUI suggests empty, wedge and symmetry from the face normals and lists the
  others for the user; pointing at a patch highlights it on the drawing.
- The material fills a box (open sides by default); the mass reported is for the mesh
  as given (its depth or wedge, wedge scaled to the full circle).
- Supported: 2D (one cell thick, empty sides) and wedge meshes about the y axis, with
  quad or triangle cross-sections, because the fill and the checker work on the
  cross-section. 3D user meshes are refused with a message (later work).

## Consequences

- `tests/test_celsetup_usermesh.sh` (runs in CI: a committed .msh, gmshToFoam): 2D tank
  imported, conditions set, 3.0 kg placed, celFoam runs; polyMesh with scale; the Couto
  wedge imported as a user mesh finds the 2.5 deg wedge and gives the same time step;
  readable errors; YAML round trip; GUI load and preview.
- The YAML reader takes one-level inline mappings (`outlet: {type: outlet, friction:
  20}`).
- celCheck's rules (outlet cells, wall cells) come from the Couto silo; on other shapes
  they are guidance, not limits.
