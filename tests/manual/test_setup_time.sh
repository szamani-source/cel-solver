#!/bin/bash
# WP 0.4 "done when" test. Run by hand in a fresh cloud session, not in CI
# (it lives outside tests/ so tests/run_tests.sh does not pick it up).
#
# From a container without OpenFOAM: tools/setup_openfoam.sh, ./Allwmake,
# tests/run_tests.sh must all succeed within 10 minutes.
#
#   tests/manual/test_setup_time.sh               # fails if OpenFOAM is already there
#   tests/manual/test_setup_time.sh --remove      # deletes /usr/lib/openfoam first
set -u
root="$(cd "$(dirname "$0")/../.." && pwd)"
limit=600

if [ "${1:-}" = --remove ]; then
    echo "removing existing OpenFOAM install and user build"
    rm -rf /usr/lib/openfoam /usr/bin/openfoam2512 "$HOME/OpenFOAM"
    find "$root" -path "*/Make/linux*" -prune -exec rm -rf {} +
fi
if [ -e /usr/lib/openfoam/openfoam2512 ]; then
    echo "OpenFOAM already installed; run with --remove to start clean"
    exit 2
fi

t0=$SECONDS
step()
{
    local name="$1"; shift
    local s=$SECONDS
    if "$@"; then
        echo "STEP  $name  $((SECONDS - s)) s"
    else
        echo "FAIL  $name after $((SECONDS - s)) s"
        exit 1
    fi
}

step setup bash "$root/tools/setup_openfoam.sh"
# OpenFOAM's bashrc is not set -u clean.
set +u
# shellcheck disable=SC1091
source /usr/lib/openfoam/openfoam2512/etc/bashrc
set -u
step build "$root/Allwmake"
step tests "$root/tests/run_tests.sh"

total=$((SECONDS - t0))
echo "TOTAL $total s (limit $limit s)"
[ "$total" -lt "$limit" ]
