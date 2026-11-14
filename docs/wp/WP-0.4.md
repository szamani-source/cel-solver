# WP 0.4: Fast OpenFOAM setup for fresh cloud sessions

State: done
Branch: `claude/project-thread-mbzxdh` (thread sessions may only push to their assigned branch)

## Goal

A new Claude cloud session can build `celFoam` and run all tests within 10 minutes of
starting, with one documented command. Today that takes about an hour, because OpenFOAM
v2512 must be compiled from source.

## Read first

- `docs/decisions/0002-openfoam-version.md`, `.github/workflows/ci.yml`
- What is known about the cloud container (found on 2026-11-11 and 2026-11-14):
  - Ubuntu 24.04, root, 4 cores. `dl.openfoam.com` is blocked by the network proxy.
  - Docker Hub image pulls returned 429 (rate limit) from inside the container.
  - The SourceForge source tarball downloads fine
    (`https://sourceforge.net/projects/openfoam/files/v2512/OpenFOAM-v2512.tgz/download`);
    a full build with `./Allwmake -j 4 -s -q -k` takes about an hour.
  - Ubuntu's own `openfoam` package is version 1912, too old.
  - GitHub Actions artifacts could not be downloaded from the container (blob host blocked).
  - A session container is reclaimed after inactivity; `/opt` does not persist.

## Build

- `tools/setup_openfoam.sh`: tries the options in order of speed and prints which worked.
  Candidates to test and time, fastest first:
  1. A prebuilt binary tarball of only what we link against (wmake, OpenFOAM, finiteVolume,
     meshTools, geometricVoF and their dependencies, blockMesh, decomposePar), published as
     a GitHub Release asset on this repository and fetched with `gh release download`.
     Check that release downloads are reachable through the proxy before building on this.
  2. A minimal from-source build of the same subset.
  3. The full from-source build (the current fallback).
- A short section in `CLAUDE.md` ("Setting up OpenFOAM in a new session") with the command.
- A decision note recording what was tried, the timings, and the choice.
- If the cloud environment supports a setup script that runs when a session starts, say so
  in the log and tell Sia how to turn it on; do not change environment settings yourself.

## Done when

- `tests/manual/test_setup_time.sh` (run by hand, not in CI) starts from a container without
  OpenFOAM, runs `tools/setup_openfoam.sh`, `./Allwmake` and `tests/run_tests.sh`, and
  finishes in under 10 minutes with all tests passing. Record the measured time in the log.
- CI is unchanged and still green.

## Out of scope

- Changing the pinned OpenFOAM version. Any solver code.

## Log

### 2026-11-14

Done:
- Probed the proxy. `mirror.gcr.io` (Google's Docker Hub mirror) serves
  `opencfd/openfoam-dev:2512` anonymously; Docker Hub with an anonymous token also worked
  today. Table of what is reachable in decision 0007.
- `tools/setup_openfoam.sh`: fetches the two OpenFOAM layers of the pinned CI image with
  curl, checks sha256, unpacks `/usr/lib/openfoam/openfoam2512` and
  `/usr/bin/openfoam2512`, apt-installs the runtime libraries. Falls back to Docker Hub,
  then to the full SourceForge build.
- Done-when test: `tests/manual/test_setup_time.sh --remove`, after also removing the
  Open MPI and Scotch packages and clearing the apt cache: setup 49 s, build 7 s,
  tests 10 s, **total 66 s**, all three tests pass. Docker Hub path alone: setup 36 s.
- The test sits in `tests/manual/`, not `tests/`, because `tests/run_tests.sh` runs every
  `tests/test_*.sh` and this one must not run in CI (it deletes OpenFOAM and calls
  `run_tests.sh` itself).
- CI workflow unchanged.

Not done, on purpose:
- GitHub Release tarball and minimal source build (reasons in decision 0007).
- The SourceForge fallback path was not timed in this package; the tarball URL still
  answers 200.

Cloud environment setup script: the environment has a "Setup script" that runs when a
session starts (Project settings, Cloud environment, gear beside the environment). Not
changed here; Sia can paste this line there so every session starts with OpenFOAM:

    curl -fsSL https://raw.githubusercontent.com/szamani-source/cel-solver/main/tools/setup_openfoam.sh | bash

Next: nothing for this package. If the mirror ever stops answering, see decision 0007.
