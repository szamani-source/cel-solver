#!/bin/bash
# WP 5.3b (decision 0073): a soft deformable footing (E 20 MPa) pushed into
# T0 shortens as the elastic column (within 20 %) and leaves the material
# density within 2 %. Builds celVoidCheck. About 1.5 min.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
bash "$root/tests/test_deformable_push.sh" --soft
