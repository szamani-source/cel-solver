# CEL solver (celFoam)

An open-source Coupled Eulerian-Lagrangian solver for large-deformation solid mechanics
(soil penetration, impact, bulk material flow), in the spirit of Abaqus/Explicit CEL and
LS-DYNA multi-material ALE.

- Eulerian side: explicit, cell-centred finite volume, built on OpenFOAM v2512 (openfoam.com).
- Lagrangian side: rigid bodies, then explicit hex8 finite elements.
- Coupling: penalty contact between Lagrangian surfaces and Eulerian material.

Status: Phase 0 (project setup). See [docs/STATUS.md](docs/STATUS.md).

## Build and test

```bash
source /path/to/OpenFOAM-v2512/etc/bashrc
./Allwmake
tests/run_tests.sh
```

## Documentation

- [docs/work-packages.md](docs/work-packages.md): the plan, split into 28 small packages
- [docs/decisions/](docs/decisions/): why things are the way they are
- [docs/spec/target-problems.md](docs/spec/target-problems.md): what the solver must solve
- [docs/reuse-survey.md](docs/reuse-survey.md): existing open-source codes and what we reuse

## Licence

GPL-3.0, because the code links against OpenFOAM (GPL-3.0). See decision 0003.
