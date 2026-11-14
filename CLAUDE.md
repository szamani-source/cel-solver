# CLAUDE.md

Read this file, then `docs/STATUS.md`, at the start of every session. Together they say
what this code is, where the work stands and what to do next. Do not rely on memory of
earlier conversations.

## What this is

`celFoam` is an explicit Coupled Eulerian-Lagrangian (CEL) solver for large-deformation
solids, like Abaqus/Explicit CEL. The plan, with reasoning, is in the project's plan doc
(link in `docs/STATUS.md`) and is broken into work packages in `docs/work-packages.md`.

## Architecture (decision 0001)

- OpenFOAM v2512 (openfoam.com) is the host only: mesh, fields, MPI decomposition,
  I/O, isoAdvector geometric VOF. **Never patch OpenFOAM's own source.**
- Our physics lives in our own libraries under `src/` with small interfaces:
  - `celConstitutive`: material models. Each model maps a strain increment and history
    to a stress at one point. Selected at runtime. Shared by Eulerian and Lagrangian sides.
  - `celEulerian`: explicit cell-centred FV solid: Lagrange step plus remap,
    N materials plus void, per-material history variables.
  - `celLagrangian`: rigid bodies first, then explicit hex8 FE with hourglass control.
  - `celCoupling`: contact search near each body, penalty forces both ways.
- `applications/solvers/celFoam` is a thin driver that calls the libraries.

One explicit cycle: stable dt from CFL, stress update, contact forces, central-difference
velocity and position update, remap onto the fixed Eulerian mesh, write.

## How to build and test

```bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc
./Allwmake
tests/run_tests.sh
```

CI (`.github/workflows/ci.yml`) runs the same in the `opencfd/openfoam-dev:2512` image.

### Setting up OpenFOAM in a new session

A fresh cloud container has no OpenFOAM. One command installs the same build CI uses,
at the same path, in about a minute (decision 0007):

```bash
tools/setup_openfoam.sh    # safe to re-run; exits at once if already installed
```

`tests/manual/test_setup_time.sh --remove` times the whole setup, build and test run
from scratch (66 s on 2026-11-14; limit 10 minutes).

## How we work (Sia's priority: small, findable pieces)

The full routine is in `docs/MANAGEMENT.md`. Follow its session start and end steps.
Bulk-solid material models (Mohr-Coulomb, Drucker-Prager) are core scope: see
`docs/constitutive-models.md`.

1. One work package per session or project thread. Pick the next one from `docs/STATUS.md`.
2. Each package has one "done when" test in `docs/work-packages.md`. Write that test first,
   as `tests/test_<name>.sh`, then make it pass. Never weaken an existing test to get green.
3. Any choice that someone might later question gets a short note in `docs/decisions/`
   (`NNNN-short-title.md`: context, decision, consequences).
4. Before finishing: update `docs/STATUS.md` (done, in progress, next, known issues),
   then commit and open a pull request named `WP x.y: <title>`.
5. Keep each package small. If it grows past about three sessions, split it and record
   the split in `docs/work-packages.md`.

## Conventions

- SI units throughout (m, kg, s, Pa). Stress positive in tension.
- OpenFOAM code style (4 spaces, Allman braces, `Info<<` for output).
- New dictionaries and fields are documented in the tutorial that first uses them.
- Tutorials live in `tutorials/`, each with `Allrun` and `Allclean`.
