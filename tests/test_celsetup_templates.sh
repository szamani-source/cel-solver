#!/bin/bash
# WP G.5: flatBottom (plane tank, axisymmetric silo) and the plane hopper
# write, mesh, fill and run 50 steps; with Gmsh, the same meshes (about 1 min).
set -e
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
python3 "$root/tests/celsetup/templates_checks.py" "$work"
