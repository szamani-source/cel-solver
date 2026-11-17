#!/bin/bash
# WP 4.0: cut-cell moving-wall Riemann coupling against penalty coupling in the 2D
# prototype (prototype/contact2d.py, decision 0042). About 30 s.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
python3 "$root/tests/contact/check_schemes.py"
