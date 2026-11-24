#!/bin/bash
# WP G.6: optional faster-time-step mesh options, off by default (about 30 s).
set -e
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
python3 "$root/tests/celsetup/fastdt_checks.py" "$work"
