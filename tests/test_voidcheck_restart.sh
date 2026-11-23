#!/bin/bash
# Decision 0116: the test driver celVoidCheck writes checkpoints
# (voidCheckDict checkpointEvery) and restarts from the latest one
# (controlDict startFrom latestTime). The coarse Couto silo (STL shell,
# nonlocal fluidity, outlet opening at 3 ms) on 2 ranks, restarted at 4 ms,
# matches the straight run at 6 ms: cells and history. About 2 minutes.
set -eu
root="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "FAIL: OpenFOAM environment not loaded"; exit 1; }

wmake "$root/tests/void/voidCheck" > "$work/log.wmake" 2>&1 \
    || { tail -n 40 "$work/log.wmake"; exit 1; }

python3 "$root/tests/restart/voidcheck_restart.py" "$work/run"
