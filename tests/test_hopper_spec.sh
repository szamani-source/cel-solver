#!/bin/bash
# WP H.1: the conical hopper validation case spec is self-consistent
# (docs/spec/conical-hopper.md, decision 0100). No OpenFOAM needed; under a second.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
python3 "$root/tests/hopper/spec_checks.py"
