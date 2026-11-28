# Installing celFoam 1.0 on Linux

celFoam is a set of libraries and programs built on top of **OpenFOAM v2512 from
openfoam.com**. Installing it means three things: get OpenFOAM v2512, load its
environment, and build celFoam against it. This page does that by hand on any Linux
distribution. Two shortcuts follow the same steps for you:

| Route | Works on | What you do |
| --- | --- | --- |
| A. By hand (this page) | any Linux, x86-64 | steps 1 to 5 below |
| B. `./install.sh` | any Linux that has OpenFOAM v2512 (step 1); with `--get-openfoam` it also installs OpenFOAM on Ubuntu 24.04 and its derivatives | `./install.sh` |
| C. Container | any Linux with Docker or Podman | `./install.sh --container` |

All three give the same `celfoam` command (step 4). Tested in CI on every change: route B
on a clean Ubuntu 24.04, as root and as a normal user, and route C with Docker
(`.github/workflows/install.yml`).

## What you need

| Item | Version | Needed for | Notes |
| --- | --- | --- | --- |
| Linux | x86-64 | everything | Ubuntu 24.04 is what CI uses. Any distribution that can run OpenFOAM v2512 works. |
| **OpenFOAM** | **v2512 (openfoam.com)** | building and running | Exactly this release. OpenFOAM 11, 12, 13 from **openfoam.org** are a different code line and will not build celFoam. Both can be installed side by side (step 2). |
| C++ compiler | g++ 9 or newer (C++17) | building | The one OpenFOAM was built with. |
| make | any | building | |
| MPI | Open MPI (the one OpenFOAM uses) | parallel runs | Serial runs need no MPI. |
| Python 3 | 3.8 or newer | tutorial scripts, tests | |
| numpy | any recent | tutorial scripts, tests | **required** |
| matplotlib | any recent | plots (`plot.py`, mesh pictures) | optional |
| ParaView | any that reads OpenFOAM | viewing results | optional |
| git | any | getting the source | optional (a tarball works too) |
| Disk, memory | about 3 GB disk, 4 GB memory | | OpenFOAM is about 2 GB. |

## Step 1. Install OpenFOAM v2512 and the system packages

Pick the section for your distribution. If there is no ready-made OpenFOAM v2512
package for it, use **From source** (works everywhere, about an hour) or route C.

### Ubuntu and Debian

OpenCFD's packages (Ubuntu and Debian releases they support):

```bash
curl -s https://dl.openfoam.com/add-debian-repo.sh | sudo bash
sudo apt-get update
sudo apt-get install openfoam2512-dev
sudo apt-get install g++ make python3 python3-numpy python3-matplotlib
```

OpenFOAM goes to `/usr/lib/openfoam/openfoam2512`.

**Ubuntu 24.04, Linux Mint 22 and other 24.04 derivatives**, shortcut: the repository
script installs the exact OpenFOAM build CI uses (pinned, checksummed, about a minute),
plus g++, make, Open MPI, python3 and numpy:

```bash
sudo tools/setup_openfoam.sh
```

### openSUSE Leap

OpenCFD's packages come from the openSUSE "science" repository:

```bash
curl -s https://dl.openfoam.com/add-science-repo.sh | sudo bash
sudo zypper refresh
sudo zypper install openfoam2512-devel
sudo zypper install gcc-c++ make python3 python3-numpy python3-matplotlib
```

If `openfoam2512-devel` is not offered for your Leap release yet, build from source.

### Fedora, RHEL, Rocky, AlmaLinux, CentOS Stream

No OpenCFD package is maintained for these at the moment; build OpenFOAM from source
(below) or use route C. The system packages for the source build:

```bash
sudo dnf install gcc-c++ make cmake flex bison zlib-devel boost-devel fftw-devel \
    readline-devel openmpi-devel python3 python3-numpy python3-matplotlib
module load mpi/openmpi-x86_64      # puts mpicc/mpirun on the PATH (add to ~/.bashrc)
```

(On RHEL-type systems some of these come from EPEL / CRB: `sudo dnf install epel-release`
and `sudo dnf config-manager --set-enabled crb`.)

### Arch Linux and Manjaro

The AUR package `openfoam-com` builds the openfoam.com release; check that its version
is v2512 before using it. Otherwise build from source. System packages:

```bash
sudo pacman -S base-devel cmake flex bison zlib boost fftw openmpi python python-numpy python-matplotlib
```

### From source (any distribution)

Install your distribution's packages for the source build (listed above for each family;
on Debian and Ubuntu: `sudo apt-get install build-essential cmake flex bison zlib1g-dev
libboost-system-dev libfftw3-dev libreadline-dev libopenmpi-dev openmpi-bin
libscotch-dev libptscotch-dev python3-numpy`), then:

```bash
mkdir -p ~/OpenFOAM && cd ~/OpenFOAM
curl -L -o OpenFOAM-v2512.tgz https://sourceforge.net/projects/openfoam/files/v2512/OpenFOAM-v2512.tgz/download
tar -xzf OpenFOAM-v2512.tgz
cd OpenFOAM-v2512
source etc/bashrc          # loads the environment of this tree
foamSystemCheck            # checks the system is ready
./Allwmake -j -s -q -l     # about an hour on 4 cores; log in log.linux64GccDPInt32Opt
```

OpenFOAM is then in `~/OpenFOAM/OpenFOAM-v2512`, and every `source .../etc/bashrc` below
uses that path. No root is needed. If `Allwmake` stops, the end of its log names the
missing package. The SourceForge page for v2512 also has a `ThirdParty-v2512.tgz`; it is
only needed when the system has no MPI or no Scotch.

## Step 2. Load the OpenFOAM environment

Every new terminal needs OpenFOAM's environment before building or running:

```bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc      # packages and tools/setup_openfoam.sh
# or
source ~/OpenFOAM/OpenFOAM-v2512/etc/bashrc            # from source
```

Check it: `echo $WM_PROJECT_VERSION` prints `v2512`, and `which blockMesh` points into
that tree.

**Another OpenFOAM on the same machine** (for example OpenFOAM 13 in `/opt/openfoam13`):
load only one per terminal. If `~/.bashrc` loads the other one, open a fresh terminal
with that line commented out, or use the `celfoam` command (step 4), which loads v2512
for you.

## Step 3. Get and build celFoam

```bash
git clone https://github.com/szamani-source/cel-solver.git   # or unpack celfoam-1.1.0.tar.gz
cd cel-solver
```

**Route B** (recommended once OpenFOAM is there): one command does the build, the
environment file, the `celfoam` command and a smoke test:

```bash
./install.sh                       # as root: /opt/celfoam-1.1.0, celfoam in /usr/local/bin
                                   # as a user: ~/celfoam-1.1.0, celfoam in ~/.local/bin
./install.sh --openfoam DIR        # OpenFOAM v2512 somewhere it does not look
./install.sh --prefix DIR          # install somewhere else
./install.sh --help
```

It finds OpenFOAM v2512 in the loaded environment or the usual places
(`/usr/lib/openfoam/openfoam2512`, `/usr/lib64/openfoam/openfoam2512`, `/opt/openfoam2512`,
`~/OpenFOAM/OpenFOAM-v2512`), checks g++, make, python3 and numpy, copies the source to
the prefix and builds there (2 to 4 minutes).

**By hand** (developers; builds in the checkout, programs go to OpenFOAM's user folders
under `~/OpenFOAM`, which are already on the path once OpenFOAM is loaded):

```bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc
./Allwmake                         # libraries, celFoam, celBodyCells
wmake tests/void/voidCheck         # celVoidCheck: the driver for the coutoHopper case
celFoam -help
```

## Step 4. Load the celFoam environment

After route B or C:

```bash
celfoam                            # a shell with OpenFOAM v2512 and celFoam loaded
celfoam <command>                  # one command, e.g. celfoam ./Allrun
source /opt/celfoam-1.1.0/etc/celfoam.sh    # the same in the current shell (your prefix)
celfoam test                       # smoke test, under a second
celfoam version
celfoam manual                     # path of the HTML manual (opens it on a desktop)
celfoam setup --gui                # the case setup page in your browser (manual: Setup tool)
celfoam check myCase               # mesh and time-step checks of a case
```

`etc/celfoam.sh` sources OpenFOAM's `etc/bashrc` and puts `<prefix>/platforms/bin` and
`<prefix>/platforms/lib` first on `PATH` and `LD_LIBRARY_PATH`. It also sets
`CELFOAM_DIR`, which the tutorials use to find their scripts.

After the by-hand build, loading OpenFOAM (step 2) is all that is needed.

## Step 5. Run a tutorial

```bash
celfoam tutorials ~/celfoam-cases          # copies all tutorials there
cd ~/celfoam-cases/solid/siloDischarge
celfoam ./Allrun                           # tank: settle, open, discharge (about 11 min)
celfoam python3 plot.py                    # snapshots.png, mass.png (needs matplotlib)
```

The tutorials are grouped by category: `solid/` (solid only), `solidFluid/` (with gas)
and `solidFluidThermal/` (with gas and temperature). The solid ones
(`tutorials/solid/README.md` lists them with run times):

| Tutorial | What | Time |
| --- | --- | --- |
| `smokeBox` | 10 x 10 x 10 elastic box, checks the install | seconds |
| `siloDischarge` | 2D tank discharging through a floor opening (T0) | about 11 min (`./Allrun 0.1`: 40 s) |
| `siloCohesion`, `siloFed` | cohesive material; a tank fed from the top while it discharges | 7 to 60 min |
| `stlSilo`, `stlHopper` | walls from STL files | 2 to 10 min |
| `coutoHopper` | the validated Couto wheat silo (axisymmetric, MC or DP) | 2 h (standard) or 8 h (fine) on 4 cores |
| `setupPlaneTank` | a plane tank written by `celfoam setup` from a short YAML file | about 1 h (`./Allrun 0.1`: minutes) |
| `stripFooting` | a rigid footing pushed into clay (coupled Eulerian-Lagrangian) | about 1 min |
| `rigidBodyMotion` | prescribed rigid bodies from STL files | seconds |

Every tutorial folder has a `README.md`. Results are standard OpenFOAM time folders: open
them in ParaView with `touch case.foam && paraview case.foam`.

Parallel runs: `celfoam mpirun -np 4 celFoam -parallel` after `decomposePar`, as in any
OpenFOAM case. As root, Open MPI needs `OMPI_ALLOW_RUN_AS_ROOT=1
OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1`; with more ranks than cores add `--oversubscribe`.

## Route C: the container

Any Linux with Docker or Podman (no OpenFOAM on the host):

```bash
./install.sh --container           # builds image celfoam:1.1.0, installs a celfoam command
celfoam tutorials cases            # same commands as above
cd cases/solid/coutoHopper && celfoam ./Allrun mc mesh
```

The image is built from the same pinned OpenFOAM v2512 image CI uses (`Dockerfile`), so
the first build downloads about 1 GB. The `celfoam` command runs each command in a fresh
container with **the current folder mounted** as its working folder, under your own user
id, so results land in your folder. Paths outside the current folder are not visible
inside; run from the case folder or a parent of it. By hand:

```bash
docker build -t celfoam:1.1.0 .
docker run --rm -it -u "$(id -u):$(id -g)" -e HOME=/tmp -v "$PWD":/work celfoam:1.1.0
```

Apptainer / Singularity (HPC clusters): `apptainer build celfoam.sif docker-daemon://celfoam:1.1.0`
on a machine with Docker, then `apptainer exec celfoam.sif /opt/celfoam/bin/celfoam ./Allrun`.

## Troubleshooting

| Message | Cause and fix |
| --- | --- |
| `OpenFOAM v2512 (openfoam.com) not found` (install.sh) | Install it (step 1), or pass `--openfoam DIR` (the folder holding `etc/bashrc`). |
| `wmake: command not found`, `WM_PROJECT_DIR: parameter null or not set` | OpenFOAM is not loaded in this terminal (step 2). |
| Build errors about unknown OpenFOAM classes | The loaded OpenFOAM is not v2512 (`echo $WM_PROJECT_VERSION`), often openfoam.org's. |
| `celFoam: command not found` | Load the environment (`celfoam`, or step 2 after a by-hand build). |
| `celVoidCheck not found` (coutoHopper) | By-hand builds: `wmake tests/void/voidCheck`. Route B builds it. |
| `No module named 'numpy'` | Install python3-numpy (step 1 for your distribution). |
| MPI error about slots, or about running as root | See Parallel runs in step 5. |
| `--get-openfoam ... needs Ubuntu 24.04` | That shortcut installs the Ubuntu 24.04 build only; elsewhere use step 1. |

The HTML manual (`docs/manual/index.html`, also `celfoam manual`) has the case setup
reference, meshing, theory and validation. Release notes: `docs/release-notes/v1.0.0.md`.
