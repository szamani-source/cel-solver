"""Runs the strip footing of tests/footing/footing_checks.py (WP 4.4) with
cell snapshots for the tutorial: python3 run_footing.py WORKDIR [--fine]."""
import os
import sys

root = os.environ.get("CEL_ROOT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
sys.path.insert(0, os.path.join(root, "tests", "footing"))
import footing_checks as fc  # noqa: E402

SNAPS = (0.2, 0.4, 0.775)      # s/B about 0.03, 0.08, 0.18
_make = fc.vc.make_case


def make_case(path, nx, ny, material, check, **kw):
    check = dict(check, writeTimes=list(SNAPS))
    return _make(path, nx, ny, material, check, **kw)


fc.vc.make_case = make_case
sys.argv = [sys.argv[0]] + sys.argv[1:]
fc.main()
