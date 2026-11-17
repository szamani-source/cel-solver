#!/bin/bash
# WP 3.6b (decision 0039), CI part: an inclined layer under gravity not along
# a mesh axis stays at rest below the stopping angle (the face reconstruction
# fix), and the mu(I) layers at E = 0.1 and 10 MPa, 5 mm, with flowImpedance
# yes, follow the Bagnold profile within 5 % with mass kept (WP 3.6c2,
# decision 0044; 0.1 MPa also against the WP 3.6c guard). About 50 s.
#
# The full WP 3.6c2 done-when (also 1 MPa; 2.5 mm and a smooth lid recorded)
# is tests/manual/test_inclined_layer.sh (about 2.5 minutes).
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/layer/inclined_layer.py" "$work" --ci
