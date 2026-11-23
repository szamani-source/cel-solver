# WP S.1 to S.7: case generator, mesh checker, setup GUI

State: in progress
Branch: `claude/project-thread-c3f4eh` (thread "Build case generator and mesh checker")

## Goal

Setting up a case takes a short YAML file (or a GUI form) instead of hand-made
blockMesh dictionaries: `celSetup case.yaml` writes the wall-aligned mesh with the
proven rules, every dictionary, the material fill and an Allrun that resumes; `celCheck`
reports the mesh against the cell-size rules before a run. Decision 0119.

## Read first

- `docs/decisions/0119-case-setup-tools.md`
- `tests/hopper/buffered_case.py` (aligned_buffered_mesh: the geometry the conical
  template reproduces), decisions 0105, 0107, 0115
- `docs/manual/meshing.html` (the cell-size rules)

## Build

- `tools/celSetup` (command), `tools/celsetup/` (package: `miniyaml.py`, `foam.py`,
  `schema.py`, `case.py`, one module per template), examples in `tools/celsetup/examples/`.

## Done when

- One test per package; see the S table in `docs/work-packages.md`.

## Out of scope

- Long runs (Sia: accuracy first, only smoke tests here). New physics.

## Log

- 2026-11-23: G.5 done. `flat.py` (blockMesh + own .geo writer for the two blocks),
  conical `shape: plane` (block_mesh_dict depth path, extruded STL, plane volume x2
  per metre). Note: the axisymmetric flat silo's step is set by the axis column (1.86
  mm effective at r 2.8 mm), the G.6 subject. celCheck's outlet rule (about 40 cells
  across) is the Couto rule and reads FAIL on small outlets; worth revisiting with Sia.
- 2026-11-23: G.4 done. `/api/open` (path or example; a relative mesh file becomes
  absolute), examples in init, applyConfig shared by start and Open, beforeunload
  guard, boundary errors marked. Next: G.5 templates.
- 2026-11-23: G.3 done (decision 0121). `usermesh.py` (import, suggestions, BTYPES,
  build, finish_import), case.py inlets and per-outlet openTime, miniyaml inline
  mappings, schema boundaries section in YAML and --describe, GUI Load mesh and
  boundaries panel (fieldset min-width 0 or the form scrolls sideways). Test fixture
  `tests/celsetup/data/tank2d.msh` (from examples/tank2d.geo) so CI tests the import
  without Gmsh. Next: G.4 polish.
- 2026-11-23: G.2 done. GUI toolbar Mesher select and "compare with the other":
  `/api/init` lists the meshers (ok + note), `/api/compare` meshes the case with the
  other mesher in its own scratch folder; overlay lines only when cells are 4 px or
  more on screen (else they hide the colours). Screenshots
  /mnt/project-files/setup/G2-gui-compare-*.png. Next: G.3 user meshes.
- 2026-11-23: G.1 done (decision 0120). Sia resumed the work as the G series (GUI,
  Gmsh or blockMesh, user meshes). `gmsh_mesh.py`: .geo from `conical.layout()`,
  Gmsh (Python module or command), gmshToFoam, patch types. Revolving drops the axis
  side from Gmsh's extrusion list (prisms on the axis); indices handle it. Couto mesh
  matches blockMesh in cells, patches and time step. Also in this step: templates
  return generic walls/outlets/bodies, plane variant of the block layout, renumbered
  our decision to 0119 (main has 0117, v1.0.0 takes 0118). Next: G.0 merge after v1.0.0, then G.2.
- 2026-11-23: S.4 done. `gui.py` + `gui/` (index.html, app.js, style.css). Endpoints
  init, geometry (instant outline), mesh (scratch case: blockMesh, fill, celCheck), save,
  run, status. Vectorised the polyMesh geometry and the fill (preview 7 s -> 2-4 s).
  Playwright screenshots checked by eye; the server test covers every call. Next: S.5
  templates.
- 2026-11-23: S.3 done. `check.py`, `tools/celCheck`. Time step reproduces
  eulerianSolid::cellSize (angled V/Smax, wedge hoop length combined) exactly. Standard
  Couto mesh: dt 5.85e-6 s from the axis cell in the spout (in-plane 2.68 mm, hoop 2.68
  mm); fine preset: 3.748e-6 s from 1.21 mm columns at the outlet edge (the shared
  columns squeeze there), not the axis. So S.6 needs both: a wider axis column and
  uniform outlet columns (edgeGrading on the hopper block's bottom edge). Next: S.4 GUI.
- 2026-11-23: S.2 done. `schema.py` holds every entry (key, default, unit, text,
  choices, basic/advanced); templates give geometry/walls/mesh, material and run are
  shared. Mesh numbers default to "preset" so `preset: fine` works from a template file.
  `celSetup --template`, `--describe`; unknown entries are errors (material passes extra
  keywords through). Empty `key:` reads as no value. Next: S.3 celCheck.
- 2026-11-23: S.1 done. `celSetup` with the conicalHopper template; the Couto example
  gives the validated 7,584-cell mesh and STL shell bit for bit, Mesh OK (non-orth max
  34.7 deg), 1503.4 kg placed, celFoam 52 steps at 5.77e-6 s, resume from the last write
  works. CI image has numpy only: own YAML reader, no matplotlib in the generator.
  Next: S.2 property schema (descriptions for every entry), then S.3 celCheck.
