# WP S.1 to S.7: case generator, mesh checker, setup GUI

State: in progress
Branch: `claude/project-thread-c3f4eh` (thread "Build case generator and mesh checker")

## Goal

Setting up a case takes a short YAML file (or a GUI form) instead of hand-made
blockMesh dictionaries: `celSetup case.yaml` writes the wall-aligned mesh with the
proven rules, every dictionary, the material fill and an Allrun that resumes; `celCheck`
reports the mesh against the cell-size rules before a run. Decision 0117.

## Read first

- `docs/decisions/0117-case-setup-tools.md`
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

- 2026-11-23: S.1 done. `celSetup` with the conicalHopper template; the Couto example
  gives the validated 7,584-cell mesh and STL shell bit for bit, Mesh OK (non-orth max
  34.7 deg), 1503.4 kg placed, celFoam 52 steps at 5.77e-6 s, resume from the last write
  works. CI image has numpy only: own YAML reader, no matplotlib in the generator.
  Next: S.2 property schema (descriptions for every entry), then S.3 celCheck.
