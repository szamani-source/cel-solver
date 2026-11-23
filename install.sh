#!/bin/bash
# celFoam installer: an optional shortcut for the steps in INSTALL.md.
#
#   ./install.sh                      # find OpenFOAM v2512, build, install to the prefix
#   ./install.sh --prefix ~/celfoam   # install somewhere else (no root needed)
#   ./install.sh --openfoam DIR       # use the OpenFOAM v2512 in DIR (holds etc/bashrc)
#   ./install.sh --get-openfoam       # Ubuntu 24.04 family only: install OpenFOAM v2512 first
#                                     #   (tools/setup_openfoam.sh, needs root)
#   ./install.sh --container          # any Linux: build the celfoam image with Docker or
#                                     #   Podman and install a `celfoam` command that runs it
#
# Other options: --jobs N (parallel build), --no-link (do not put `celfoam` on the PATH),
# --no-test (skip the smoke test), --help.
#
# What it does (native route, decision 0118):
#   1. finds OpenFOAM v2512 (openfoam.com): --openfoam, the loaded environment, or the
#      usual install paths; stops with a pointer to INSTALL.md if there is none
#   2. copies this source tree to the prefix (default /opt/celfoam-<version> as root,
#      ~/celfoam-<version> otherwise)
#   3. builds the libraries, celFoam, celBodyCells and celVoidCheck into <prefix>/platforms
#   4. writes <prefix>/etc/celfoam.sh (the environment) and <prefix>/bin/celfoam (the
#      command), links `celfoam` into /usr/local/bin (root) or ~/.local/bin
#   5. runs the smoke test (tests/test_smokeBox.sh, under a second)
set -euo pipefail

src="$(cd "$(dirname "$0")" && pwd)"
version="$(cat "$src/VERSION" 2>/dev/null || echo dev)"
prefix=""
foam=""
getFoam=no
container=no
link=yes
smoke=yes
jobs="$(nproc 2>/dev/null || echo 2)"

usage() { sed -n '2,24p' "$0" | sed 's/^# \{0,1\}//'; }
log() { echo "[celfoam-install] $*"; }
die() { echo "[celfoam-install] ERROR: $*" >&2; exit 1; }

while [ $# -gt 0 ]; do
    case "$1" in
        --prefix) prefix="$2"; shift 2 ;;
        --openfoam) foam="$2"; shift 2 ;;
        --get-openfoam) getFoam=yes; shift ;;
        --container) container=yes; shift ;;
        --jobs) jobs="$2"; shift 2 ;;
        --no-link) link=no; shift ;;
        --no-test) smoke=no; shift ;;
        -h|--help) usage; exit 0 ;;
        *) die "unknown option $1 (see --help)" ;;
    esac
done

SUDO=""
isRoot=no
[ "$(id -u)" -eq 0 ] && isRoot=yes

# ---------------------------------------------------------------- container route
if [ "$container" = yes ]; then
    engine=""
    for e in docker podman; do
        if command -v "$e" > /dev/null; then engine=$e; break; fi
    done
    [ -n "$engine" ] || die "--container needs docker or podman (see INSTALL.md, route C)"
    image="celfoam:$version"
    log "building image $image with $engine (downloads the OpenFOAM v2512 image, about 1 GB)"
    "$engine" build -t "$image" -t celfoam:latest "$src"
    if [ "$isRoot" = yes ]; then bindir=/usr/local/bin; else bindir="$HOME/.local/bin"; fi
    mkdir -p "$bindir"
    cat > "$bindir/celfoam" <<EOF
#!/bin/sh
# celfoam $version, container route (install.sh --container). Runs the celfoam command
# inside the image $image with the current folder mounted as /work.
tty=""; [ -t 0 ] && [ -t 1 ] && tty="-t"
exec $engine run --rm -i \$tty -u "\$(id -u):\$(id -g)" -e HOME=/tmp \\
    -v "\$PWD":/work -w /work $image "\$@"
EOF
    chmod +x "$bindir/celfoam"
    log "installed $bindir/celfoam (runs $image)"
    case ":$PATH:" in *":$bindir:"*) ;; *) log "add $bindir to your PATH" ;; esac
    if [ "$smoke" = yes ]; then
        "$bindir/celfoam" test
    fi
    exit 0
fi

# ---------------------------------------------------------------- native route
# 1. OpenFOAM v2512
isFoam() { [ -f "$1/etc/bashrc" ] && grep -qE "api=2512|WM_PROJECT_VERSION=v2512" "$1/META-INFO/api-info" "$1/etc/bashrc" 2>/dev/null; }

findFoam()
{
    local d
    for d in "${WM_PROJECT_DIR:-}" \
             /usr/lib/openfoam/openfoam2512 /usr/lib64/openfoam/openfoam2512 \
             /opt/openfoam2512 /opt/OpenFOAM/OpenFOAM-v2512 /opt/OpenFOAM-v2512 \
             "$HOME/OpenFOAM/OpenFOAM-v2512" "$HOME/openfoam/OpenFOAM-v2512"; do
        if [ -n "$d" ] && isFoam "$d"; then echo "$d"; return 0; fi
    done
    return 1
}

if [ -z "$foam" ]; then
    foam="$(findFoam || true)"
fi
if [ -z "$foam" ] && [ "$getFoam" = yes ]; then
    . /etc/os-release 2>/dev/null || true
    if [ "${UBUNTU_CODENAME:-${VERSION_CODENAME:-}}" != noble ]; then
        die "--get-openfoam installs the Ubuntu 24.04 build and needs Ubuntu 24.04 (or a
    derivative such as Mint 22). On this system install OpenFOAM v2512 as INSTALL.md
    describes for your distribution, or use --container."
    fi
    log "installing OpenFOAM v2512 with tools/setup_openfoam.sh"
    "$src/tools/setup_openfoam.sh"
    foam=/usr/lib/openfoam/openfoam2512
fi
if [ -z "$foam" ]; then
    die "OpenFOAM v2512 (openfoam.com) not found. Install it first (INSTALL.md, step 1),
    then run again, with --openfoam DIR if it is not in a standard place.
    Shortcuts: --get-openfoam (Ubuntu 24.04 family) or --container (any Linux)."
fi
isFoam "$foam" || die "$foam does not look like OpenFOAM v2512 (no etc/bashrc for v2512)"
log "OpenFOAM v2512: $foam"

# 2. prerequisites that the build and the scripts need
missing=""
command -v python3 > /dev/null || missing="$missing python3"
python3 -c 'import numpy' 2>/dev/null || missing="$missing python3-numpy"
for t in g++ make; do command -v "$t" > /dev/null || missing="$missing $t"; done
[ -z "$missing" ] || die "missing:$missing (INSTALL.md, step 1 lists the packages per distribution)"
command -v mpirun > /dev/null || log "note: mpirun not on the PATH yet; OpenFOAM's environment may add it. Parallel runs need MPI."

# 3. copy the source tree to the prefix
if [ -z "$prefix" ]; then
    if [ "$isRoot" = yes ]; then prefix="/opt/celfoam-$version"; else prefix="$HOME/celfoam-$version"; fi
fi
prefix="$(mkdir -p "$prefix" 2>/dev/null && cd "$prefix" && pwd || echo "$prefix")"
if [ ! -w "$(dirname "$prefix")" ] && [ ! -w "$prefix" ]; then
    command -v sudo > /dev/null || die "cannot write $prefix; choose another --prefix"
    SUDO=sudo
fi
log "installing celFoam $version to $prefix"
$SUDO mkdir -p "$prefix"
if [ "$src" != "$prefix" ]; then
    (cd "$src" && tar -cf - \
        --exclude=.git --exclude='*/Make/linux64*' --exclude=lnInclude --exclude='*.o' \
        --exclude='*.dep' --exclude='tests/*.log' --exclude=__pycache__ --exclude=platforms \
        --exclude='tutorials/*/log.*' --exclude='tutorials/*/processor*' --exclude='tutorials/*/run' \
        .) | $SUDO tar -xf - -C "$prefix"
fi

# 4. build into <prefix>/platforms
log "building (about 2 to 4 minutes; log in $prefix/log.build)"
$SUDO bash -c "
    set +eu
    source '$foam/etc/bashrc'
    set -e
    export FOAM_USER_APPBIN='$prefix/platforms/bin' FOAM_USER_LIBBIN='$prefix/platforms/lib'
    export WM_NCOMPPROCS='$jobs'
    export LD_LIBRARY_PATH=\"\$FOAM_USER_LIBBIN:\$LD_LIBRARY_PATH\"
    mkdir -p \"\$FOAM_USER_APPBIN\" \"\$FOAM_USER_LIBBIN\"
    cd '$prefix'
    ./Allwmake
    wmake tests/void/voidCheck
" > /tmp/celfoam-build.$$.log 2>&1 || { tail -40 /tmp/celfoam-build.$$.log; die "build failed (full log /tmp/celfoam-build.$$.log)"; }
$SUDO cp /tmp/celfoam-build.$$.log "$prefix/log.build"
rm -f /tmp/celfoam-build.$$.log
for b in celFoam celBodyCells celVoidCheck; do
    [ -x "$prefix/platforms/bin/$b" ] || die "build did not produce $b (see $prefix/log.build)"
done

# 5. environment file and command
tmp="$(mktemp -d)"
cat > "$tmp/celfoam.sh" <<EOF
# celFoam $version environment (written by install.sh). Load it with
#     source $prefix/etc/celfoam.sh
# It loads OpenFOAM v2512 and puts the celFoam programs and libraries first.
CELFOAM_DIR='$prefix'
CELFOAM_VERSION='$version'
CELFOAM_OPENFOAM='$foam'
export CELFOAM_DIR CELFOAM_VERSION CELFOAM_OPENFOAM
if [ "\${WM_PROJECT_DIR:-}" != "\$CELFOAM_OPENFOAM" ]; then
    . "\$CELFOAM_OPENFOAM/etc/bashrc" || true
fi
export PATH="\$CELFOAM_DIR/platforms/bin:\$PATH"
export LD_LIBRARY_PATH="\$CELFOAM_DIR/platforms/lib\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}"
EOF
$SUDO mkdir -p "$prefix/etc" "$prefix/bin"
$SUDO cp "$tmp/celfoam.sh" "$prefix/etc/celfoam.sh"
$SUDO cp "$src/tools/celfoam" "$prefix/bin/celfoam"
$SUDO chmod 755 "$prefix/bin/celfoam"
rm -rf "$tmp"

if [ "$link" = yes ]; then
    if [ -w /usr/local/bin ] || [ "$SUDO" = sudo ]; then
        $SUDO ln -sf "$prefix/bin/celfoam" /usr/local/bin/celfoam && bindir=/usr/local/bin
    else
        mkdir -p "$HOME/.local/bin"
        ln -sf "$prefix/bin/celfoam" "$HOME/.local/bin/celfoam"
        bindir="$HOME/.local/bin"
    fi
    log "linked $bindir/celfoam"
    case ":$PATH:" in *":$bindir:"*) ;; *) log "add $bindir to your PATH" ;; esac
fi

# 6. smoke test
if [ "$smoke" = yes ]; then
    "$prefix/bin/celfoam" test
fi

cat <<EOF

celFoam $version is installed in $prefix.

  celfoam                      a shell with OpenFOAM v2512 and celFoam loaded
  celfoam tutorials [DIR]      copy the tutorials to DIR (default ./celfoam-tutorials)
  celfoam manual               where the HTML manual is
  source $prefix/etc/celfoam.sh   the same environment in the current shell

Next: INSTALL.md, step 5 (run a tutorial).
EOF
