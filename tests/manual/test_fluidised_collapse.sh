#!/bin/bash
# WP F.8: initially fluidised dam break in air against Roche et al. (2008),
# dry collapse with gas off unchanged, 2 and 4 processors, cost (decision
# 0133). tests/gas/collapse_checks.py; about 8 minutes on 4 cores.
#
# Checkpointed: set CEL_WORK to a lasting folder (for example under
# /mnt/project-files) and re-run after an interruption; every run resumes
# from its last write (every 0.02 s).
set -eu
root="$(cd "$(dirname "$0")/../.." && pwd)"
work="${CEL_WORK:-$(mktemp -d)}"
[ -n "${CEL_WORK:-}" ] || trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

python3 "$root/tests/gas/collapse_checks.py" "$work" "$@"
