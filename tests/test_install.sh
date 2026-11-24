#!/bin/bash
# WP L.3 (decision 0118): install.sh installs a working copy. Installs this checkout to
# a temporary prefix against the loaded OpenFOAM v2512 (no link, no root), then runs
# from the installed copy: the smoke test, the tutorial copy, smokeBox's Allrun and the
# coutoHopper mesh check. The clean-system installs (fresh Ubuntu 24.04, container)
# run in .github/workflows/install.yml. About 4 minutes.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
prefix="$work/celfoam"
foam="${WM_PROJECT_DIR:-/usr/lib/openfoam/openfoam2512}"

"$root/install.sh" --prefix "$prefix" --openfoam "$foam" --no-link
cf="$prefix/bin/celfoam"
"$cf" version | grep -q "celFoam $(cat "$root/VERSION")"
"$cf" which celFoam | grep -q "^$prefix/platforms/bin/celFoam$"
"$cf" which celVoidCheck | grep -q "^$prefix/platforms/bin/"
test -f "$("$cf" manual)"

cd "$work"
"$cf" tutorials cases > /dev/null
(cd cases/smokeBox && "$cf" ./Allrun > /dev/null && grep -q '^End' log.celFoam)
(cd cases/coutoHopper && "$cf" ./Allrun mc mesh > mesh.log 2>&1 \
    && grep -q "Mesh OK" mesh.log && grep -q "cells: *7584" mesh.log)
(cd cases/setupPlaneTank && "$cf" setup planeTank.yaml case > setup.log 2>&1 \
    && { "$cf" check case > check.log 2>&1 || true; } && grep -q "PASS  time step" check.log)
"$cf" setup --describe flatBottom | grep -q outletDiameter
echo "install test: PASS"
