#!/bin/bash
# Install OpenFOAM v2512 (openfoam.com) in a fresh Ubuntu 24.04 container, fast.
#
#   tools/setup_openfoam.sh            # then: source /usr/lib/openfoam/openfoam2512/etc/bashrc
#
# Result: the same tree as the CI image opencfd/openfoam-dev:2512, at the same path
# (/usr/lib/openfoam/openfoam2512). Methods, fastest first (decision 0006):
#
#   1. registry: download the OpenFOAM layers of the pinned opencfd/openfoam-dev:2512
#      image straight from a registry (no docker needed), check their sha256, unpack only
#      /usr/lib/openfoam and /usr/bin/openfoam2512, and apt-install the Ubuntu runtime
#      libraries. About 1 minute. Registries tried: mirror.gcr.io (Google's Docker Hub
#      mirror, no rate limit seen), then registry-1.docker.io (rate limited at times).
#   2. source: full build from the SourceForge tarball. About an hour on 4 cores.
#
# Needs root (or sudo), curl, python3, tar. Safe to run again: it exits early if
# OpenFOAM is already installed. Set CEL_OF_METHOD=registry or =source to force one,
# CEL_OF_REGISTRIES to change the registry list.
set -euo pipefail

prefix=/usr/lib/openfoam
foamDir=$prefix/openfoam2512
repo=opencfd/openfoam-dev
# amd64 manifest of opencfd/openfoam-dev:2512 (image index sha256:780d7c1c...),
# pinned so every session gets the same bits. To move, see decision 0006.
manifest=${CEL_OF_MANIFEST:-sha256:0c279c90cc0531e4f5acfe5365c529d5794a50fc7d8e4a0a45037c5f404f21b9}
registries=${CEL_OF_REGISTRIES:-"mirror.gcr.io registry-1.docker.io"}
sfUrl=https://sourceforge.net/projects/openfoam/files/v2512/OpenFOAM-v2512.tgz/download

# Ubuntu packages the OpenFOAM debs depend on, minus compilers already in the image and
# minus what is only needed to rebuild OpenFOAM itself (gfortran, CGAL, boost headers).
runtimeDeps="mpi-default-bin mpi-default-dev libboost-system1.83.0 libfftw3-double3
    libgmp10 libmpfr6 libptscotch-7.0 libscotch-7.0 libreadline8t64 flex libfl-dev
    libscotch-dev libptscotch-dev zlib1g-dev g++ make python3"
# Extra packages a from-source build needs.
sourceDeps="libreadline-dev libboost-system-dev libcgal-dev libfftw3-dev gfortran"

SUDO=""
[ "$(id -u)" -eq 0 ] || SUDO=sudo

log() { echo "[setup_openfoam $(date +%H:%M:%S)] $*"; }

installed()
{
    [ -f "$foamDir/etc/bashrc" ] \
    && [ -f "$foamDir/platforms/linux64GccDPInt32Opt/lib/libOpenFOAM.so" ]
}

aptInstall()
{
    log "apt-get install $*"
    $SUDO apt-get update -qq
    DEBIAN_FRONTEND=noninteractive $SUDO apt-get install -y -qq --no-install-recommends "$@" \
        > /tmp/setup_openfoam_apt.log 2>&1 \
        || { tail -20 /tmp/setup_openfoam_apt.log; return 1; }
}

ensureNumpy()
{
    python3 -c 'import numpy' 2>/dev/null || aptInstall python3-numpy
}

# Anonymous bearer token for a registry, or nothing if it needs none.
registryToken()
{
    case "$1" in
        registry-1.docker.io)
            curl -fsS "https://auth.docker.io/token?service=registry.docker.io&scope=repository:$repo:pull" \
            | python3 -c 'import json,sys; print(json.load(sys.stdin)["token"])'
            ;;
    esac
}

fromRegistry()
{
    local reg="$1" work token auth=()
    work="$(mktemp -d)"
    token="$(registryToken "$reg")" || return 1
    [ -z "$token" ] || auth=(-H "Authorization: Bearer $token")

    log "fetching manifest $manifest from $reg"
    curl -fsSL "${auth[@]}" -H 'Accept: application/vnd.oci.image.manifest.v1+json' \
        "https://$reg/v2/$repo/manifests/$manifest" -o "$work/manifest.json" || return 1
    echo "${manifest#sha256:}  $work/manifest.json" | sha256sum -c --quiet || return 1

    # Layers of the image, in order. The OpenFOAM tree is in the two apt-install layers
    # (224 MB: openfoam2512, 178 MB: openfoam2512-dev). The Ubuntu base (30 MB) and the
    # small config layers are skipped by size; extraction below only takes OpenFOAM paths.
    local layers
    layers=$(python3 -c '
import json,sys
m=json.load(open(sys.argv[1]))
print(" ".join(l["digest"] for l in m["layers"] if l["size"] > 50*2**20))' \
        "$work/manifest.json")

    log "downloading layers: $layers"
    local d pids=()
    for d in $layers; do
        curl -fsSL "${auth[@]}" "https://$reg/v2/$repo/blobs/$d" -o "$work/${d#sha256:}.tgz" &
        pids+=($!)
    done
    for d in "${pids[@]}"; do wait "$d" || return 1; done
    for d in $layers; do
        echo "${d#sha256:}  $work/${d#sha256:}.tgz" | sha256sum -c --quiet || return 1
    done

    aptInstall $runtimeDeps

    log "unpacking OpenFOAM into $prefix"
    $SUDO mkdir -p "$prefix"
    for d in $layers; do
        # Apply layers in order; not every layer holds both paths.
        $SUDO tar -xzf "$work/${d#sha256:}.tgz" -C / --wildcards \
            'usr/lib/openfoam/*' 'usr/bin/openfoam2512' 2>/dev/null || true
    done
    rm -rf "$work"
    installed
}

fromSource()
{
    log "from-source build: this takes about an hour on 4 cores"
    aptInstall $runtimeDeps $sourceDeps
    local work
    work="$(mktemp -d)"
    curl -fsSL "$sfUrl" -o "$work/OpenFOAM-v2512.tgz"
    $SUDO mkdir -p "$prefix"
    $SUDO tar -xzf "$work/OpenFOAM-v2512.tgz" -C "$prefix"
    $SUDO rm -rf "$foamDir"
    $SUDO mv "$prefix/OpenFOAM-v2512" "$foamDir"
    rm -rf "$work"
    (
        cd "$foamDir"
        set +eu
        # shellcheck disable=SC1091
        source etc/bashrc WM_MPLIB=SYSTEMOPENMPI
        ./Allwmake -j "$(nproc)" -s -q -k > "$foamDir/log.Allwmake" 2>&1
    )
    installed
}

main()
{
    local t0=$SECONDS
    if installed; then
        log "OpenFOAM v2512 already at $foamDir; nothing to do"
        ensureNumpy
        return 0
    fi

    local method=${CEL_OF_METHOD:-auto} ok=""
    if [ "$method" = auto ] || [ "$method" = registry ]; then
        for reg in $registries; do
            if fromRegistry "$reg"; then ok="registry ($reg)"; break; fi
            log "registry $reg failed, trying next"
        done
    fi
    if [ -z "$ok" ] && { [ "$method" = auto ] || [ "$method" = source ]; }; then
        fromSource && ok="source (SourceForge tarball)"
    fi
    if [ -z "$ok" ]; then
        log "FAILED: no method worked"
        return 1
    fi

    ensureNumpy
    log "done in $((SECONDS - t0)) s via $ok"
    log "next: source $foamDir/etc/bashrc"
}

main "$@"
