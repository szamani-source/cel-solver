#!/bin/bash
# WP G.3: celSetup takes a user-provided mesh (.msh, polyMesh, a wedge),
# assigns boundary conditions to its patches, fills and runs it (about 40 s).
set -e
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
python3 "$root/tests/celsetup/usermesh_checks.py" "$work"
