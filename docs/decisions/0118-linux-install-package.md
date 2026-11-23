# 0118: Version 1.0 and how it installs on Linux

Date: 2026-11-23. Thread "Package version 1 for Linux installation" (WP L.1-L.5).

## Context

Sia (2026-11-23): "We are ready to deploy the first version for linux installation", it
"should be ok for any linux", with a read-me that says what is needed and how to source
OpenFOAM. celFoam builds against OpenFOAM v2512 (openfoam.com) only (decision 0002), and
OpenFOAM v2512 is not packaged the same way everywhere: OpenCFD ships Debian/Ubuntu
packages and openSUSE packages, no maintained Fedora/RHEL package, an AUR build for Arch,
and the source tarball. The laptop install (2026-11-21) used `tools/setup_openfoam.sh`,
which unpacks the Ubuntu 24.04 build of the CI image.

Nonlocal mesh dependence, cohesion and compaction were listed by Sia as items before the
Linux install; they are in the code (decisions 0111-0114) with open gaps. Sia confirmed
2026-11-23 that they ship as known limitations of version 1. The setup tool (PR #77) is
not part of it.

## Decision

1. **The README is the primary route.** `INSTALL.md` gives the steps by hand for any
   Linux: OpenFOAM v2512 per distribution family (OpenCFD packages where they exist,
   otherwise from source, with the system packages for each), sourcing
   `etc/bashrc`, building, loading the environment, running a tutorial. The manual's
   installation page summarises it.
2. **`install.sh` is a convenience that follows the same steps.** It does not install
   OpenFOAM by itself (that is distribution specific) except with `--get-openfoam` on
   the Ubuntu 24.04 family, where it runs `tools/setup_openfoam.sh`. It finds an
   existing v2512, copies the source to a prefix (`/opt/celfoam-<version>` as root,
   `~/celfoam-<version>` as a user), builds there with `FOAM_USER_APPBIN/LIBBIN` set to
   `<prefix>/platforms`, writes `<prefix>/etc/celfoam.sh` and the `celfoam` command
   (`tools/celfoam`), and runs the smoke test. Copying the whole tree (tests included)
   keeps the tutorials' scripts and the manual next to the programs, and lets the
   coutoHopper tutorial use `tests/hopper/` through `$CELFOAM_DIR`.
3. **The container covers every other system.** `Dockerfile` starts from the same
   pinned `opencfd/openfoam-dev:2512` manifest as CI and `setup_openfoam.sh` (decision
   0006) and runs `install.sh` inside. `install.sh --container` builds it with Docker or
   Podman and installs a `celfoam` wrapper that mounts the current folder, so the
   commands are the same as natively. The image is built locally, not published to a
   registry, to keep the release to one artefact.
4. **No .deb or .rpm.** A package would still need OpenCFD's repository for its
   OpenFOAM dependency, would be tied to one distribution's library versions, and
   would need one package per family. The source tarball plus a build of 2 to 4 minutes
   works everywhere OpenFOAM does.
5. **Versioning and release.** `VERSION` holds the version (semantic versioning: the
   major number changes when case files stop working unchanged). On main, when `VERSION`
   names a version with no tag, `.github/workflows/release.yml` tags the commit
   `v<version>` and publishes a GitHub release with `git archive` tarball, its SHA-256
   and `docs/release-notes/v<version>.md`. `CHANGELOG.md` lists changes per version.
   `LICENSE` holds the GPL-3.0 text (decision 0003).
6. **Tests.** `tests/test_install.sh` (in the regression suite): install to a temporary
   prefix and run from it. `.github/workflows/install.yml`: a clean Ubuntu 24.04
   container with only curl, python3, tar and sudo, `install.sh --get-openfoam` as root,
   then a second install as a normal user into the home folder; and
   `install.sh --container` with Docker. Native installs on other distributions are not
   tested in CI (no OpenFOAM v2512 package to install from there); the container route
   is.

## Consequences

- Users need OpenFOAM v2512 first; on distributions without a package that is an hour
  of source build, or the container.
- The installed copy is read-only for users when installed as root; `celfoam tutorials`
  copies cases out before running them.
- `tutorials/coutoHopper` depends on the scripts in `tests/hopper/`; moving them breaks
  the tutorial (the install test catches it).
- Version 2 (gas in the void, decision to come) bumps `VERSION` and adds release notes;
  the release workflow does the rest.
