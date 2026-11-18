#!/bin/bash
# WP 4.3e at 5 mm cells (pile 20 cells, about 25 min): rough P / (cu D) 11.26
# against 11.94 (-5.7 %) on 2026-11-18.
set -eu
root="$(cd "$(dirname "$0")/../.." && pwd)"
bash "$root/tests/test_rough_pile.sh" --fine
