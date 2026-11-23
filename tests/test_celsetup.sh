#!/bin/bash
# WP S.1: celSetup writes a complete, runnable Couto hopper case (about 30 s).
set -e
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
python3 "$root/tests/celsetup/celsetup_checks.py" "$work"
