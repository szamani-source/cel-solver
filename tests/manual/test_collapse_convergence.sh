#!/bin/bash
# WP F.8k done-when (decision 0138): dry dam break, plain Mohr-Coulomb with
# yieldShearJump no and slowFlowJump yes, run-out (front at 20 mm thickness)
# on 10, 5 and 2.5 mm cells within 5 % of the 2.5 mm one at Courant 0.4. About 2.5 h on 4 cores.
#
#   tests/manual/test_collapse_convergence.sh            # 10 5 2.5 mm
#   tests/manual/test_collapse_convergence.sh 10 5       # chosen meshes
#   tests/manual/test_collapse_convergence.sh --base 10  # without slowFlowJump
set -eu
root="$(cd "$(dirname "$0")/../.." && pwd)"
work="${WORK:-$(mktemp -d)}"
[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }
python3 "$root/tests/collapse/collapse_convergence.py" "$work" "$@"
