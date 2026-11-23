#!/bin/bash
# WP S.4: the celSetup GUI server answers every call the page makes (about 15 s).
set -e
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
python3 "$root/tests/celsetup/gui_checks.py" "$work"
