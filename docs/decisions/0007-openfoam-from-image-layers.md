# 0007 Install OpenFOAM in cloud sessions from the CI image's layers

Date: 2026-11-14. Status: accepted (WP 0.4).

## Context

A fresh Claude cloud session (Ubuntu 24.04, 4 cores, root, egress through a proxy) has no
OpenFOAM. Building v2512 from the SourceForge tarball takes about an hour. Known limits:
`dl.openfoam.com` (the .deb repository) is blocked, `docker pull` from Docker Hub hit 429
rate limits, and there is no Docker daemon anyway.

The brief listed three candidates: a prebuilt tarball as a GitHub Release asset, a
minimal from-source build, and the full from-source build.

What was found on 2026-11-14:

| Source | Reachable | Notes |
| --- | --- | --- |
| `mirror.gcr.io` (Google's Docker Hub mirror) | yes, anonymous, no token | serves `opencfd/openfoam-dev:2512` |
| `registry-1.docker.io` with anonymous token | yes (this time) | rate limited per IP at times |
| `archive.ubuntu.com` (apt) | yes | |
| SourceForge tarball | yes | full build about 1 h |
| `dl.openfoam.com`, `develop.openfoam.com`, `quay.io` | no | blocked |

The CI image is Ubuntu 24.04 too. Its two `apt-get install` layers (224 MB and 178 MB
compressed) hold the openfoam2512 and openfoam2512-dev packages, which put everything
under `/usr/lib/openfoam/openfoam2512` plus the wrapper `/usr/bin/openfoam2512`.
A container registry is plain HTTPS, so these layers can be fetched with `curl` and
unpacked with `tar`; no Docker is needed.

## Decision

`tools/setup_openfoam.sh` downloads those two layers of the pinned amd64 manifest
(`sha256:0c279c90...`), checks their sha256, unpacks only the OpenFOAM paths, and
apt-installs the Ubuntu libraries the openfoam2512 packages depend on (Open MPI, Scotch,
FFTW, flex and a few more). Registries in order: `mirror.gcr.io`, then Docker Hub. The
full SourceForge build stays as the last fallback.

Timings from a container with no OpenFOAM and an empty apt cache:

| Step | Time |
| --- | --- |
| download two layers (400 MB) | 10 s |
| apt-get update and install | 19 s |
| unpack | 18 s |
| `./Allwmake` (celFoam from clean) | 7 s |
| `tests/run_tests.sh` | 10 s |
| total (`tests/manual/test_setup_time.sh`) | 66 s |

The GitHub Release option was not built: it would need a tarball that someone keeps up to
date, while the image layers are maintained by OpenCFD and are bit-identical to what CI
uses. The minimal from-source build was not built either: it would still take tens of
minutes and is not needed while a registry is reachable.

## Consequences

- Cloud sessions and CI run the same OpenFOAM build, at the same path.
- The image is pinned by manifest digest. If OpenCFD pushes a new `2512` tag, CI follows it
  but the setup script does not. To move: read the new amd64 digest from
  `https://mirror.gcr.io/v2/opencfd/openfoam-dev/manifests/2512` and update `manifest=` in
  the script; `CEL_OF_MANIFEST` overrides it for a trial.
- The script picks layers by size (over 50 MB). If a new image changes its layer layout,
  check that rule.
- Depends on `mirror.gcr.io` or Docker Hub staying reachable from the sandbox. If both go,
  the script falls back to the one-hour source build on its own (not timed in WP 0.4).
- Only amd64 is handled.
