# Changelog

All notable changes to celFoam. Versions follow semantic versioning: the major number
changes when case files stop working unchanged.

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
