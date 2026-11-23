# Changelog

All notable changes to celFoam. Versions follow semantic versioning: the major number
changes when case files stop working unchanged.

## 1.1.0 (2026-11-23)

Case setup tools. Release notes: `docs/release-notes/v1.1.0.md`. Case files from 1.0.0
run unchanged.

- `celfoam setup` (celSetup): a complete case from a short YAML file or from a setup page
  in the browser (`celfoam setup --gui`); templates conicalHopper (wedge or plane),
  flatBottom (silo or plane tank) and userMesh (your own mesh, a boundary condition per
  patch); meshing with blockMesh or Gmsh (decisions 0119-0121).
- `celfoam check` (celCheck): time step as celFoam computes it, outlet cells in grain
  diameters, walls, non-orthogonality, size jumps, run time.
- Optional faster-time-step mesh options, off and not validated (decision 0122).
- Manual page "Setup tool", tutorial `setupPlaneTank`.

## 1.0.0 (2026-11-23)

First release. Release notes: `docs/release-notes/v1.0.0.md`.

- Eulerian solid with free surface, outlets and inlets; Lagrangian rigid and hex8
  bodies; penalty contact; STL walls; axisymmetric wedge; wall-aligned angled meshes;
  MPI; restart and checkpoints.
- Materials: linear elastic, von Mises, Mohr-Coulomb, Drucker-Prager, with cohesion,
  dilation, mu(I), nonlocal fluidity, compaction cap and K0 relaxation.
- Installation for any Linux (`INSTALL.md`, `install.sh`, `Dockerfile`, decision 0118),
  the `celfoam` command, install tests in CI, releases from `VERSION`.
- Tutorials, including the validated Couto silo (`tutorials/coutoHopper`).
- HTML user manual (`docs/manual/`).
