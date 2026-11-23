#!/bin/bash
# WP G.1: celSetup meshes the Couto hopper with Gmsh as well as blockMesh
# (about 40 s with Gmsh; the meshing part is skipped where Gmsh is absent).
set -e
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
python3 "$root/tests/celsetup/gmsh_checks.py" "$work"
