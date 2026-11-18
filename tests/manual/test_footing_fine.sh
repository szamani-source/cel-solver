#!/bin/bash
# WP 4.4 at 5 mm cells (B = 20 cells, about 9 min): q / cu 5.04 on 2026-11-18.
set -eu
root="$(cd "$(dirname "$0")/../.." && pwd)"
bash "$root/tests/test_footing.sh" --fine
