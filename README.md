# CEL solver (celFoam)

An open-source Coupled Eulerian-Lagrangian solver for large-deformation solid mechanics
(soil penetration, impact, bulk material flow), in the spirit of Abaqus/Explicit CEL and
LS-DYNA multi-material ALE.

- Eulerian side: explicit, cell-centred finite volume, built on OpenFOAM v2512 (openfoam.com).
- Lagrangian side: rigid bodies, then explicit hex8 finite elements.
- Coupling: penalty contact between Lagrangian surfaces and Eulerian material.

Version 1.1.0 (2026-11-23, adds the case setup tools `celfoam setup` and `celfoam check`): release notes in
[docs/release-notes/v1.1.0.md](docs/release-notes/v1.1.0.md), changes in
[CHANGELOG.md](CHANGELOG.md), work status in [docs/STATUS.md](docs/STATUS.md).

## Install

Any Linux with OpenFOAM v2512 (openfoam.com), or any Linux with Docker or Podman.
Step by step, with the packages per distribution: **[INSTALL.md](INSTALL.md)**. Short form:

```bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc   # OpenFOAM v2512 loaded
./install.sh                                       # build, install, smoke test
celfoam tutorials ~/cases && cd ~/cases/solid/siloDischarge && celfoam ./Allrun
```

or `./install.sh --container` without OpenFOAM on the host.

## Build and test (developers)

```bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc
./Allwmake
tests/run_tests.sh
```

## Documentation

- [docs/manual/index.html](docs/manual/index.html): the user manual (HTML, works offline)

- [docs/work-packages.md](docs/work-packages.md): the plan, split into 28 small packages
- [docs/decisions/](docs/decisions/): why things are the way they are
- [docs/spec/target-problems.md](docs/spec/target-problems.md): what the solver must solve
- [docs/reuse-survey.md](docs/reuse-survey.md): existing open-source codes and what we reuse

## Licence

GPL-3.0 ([LICENSE](LICENSE)), because the code links against OpenFOAM (GPL-3.0). See decision 0003.
