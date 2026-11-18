#!/bin/bash
# WP 5.3b (not passing yet, decision 0072): soft deformable footing (E 20 MPa)
# pushed into T0: shortening against the elastic column and material density.
# On 2026-11-18: ratio 1.44 (band 20 %), density 1.62 rho0 (band 1.02). About 1.5 min.
set -eu
root="$(cd "$(dirname "$0")/../.." && pwd)"
bash "$root/tests/test_deformable_push.sh" --soft
