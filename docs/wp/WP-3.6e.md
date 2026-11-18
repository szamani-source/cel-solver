# WP 3.6e: column collapse with mu(I) against experiments

State: done at 10 mm (2026-11-18)
Branch: `claude/project-thread-4zdz5u` (with WP 3.6g and 3.6h, PR #44)

## Goal

Granular column collapse with mu(I) in the decision 0019 bands at a = 0.5 and 2, at
T0's real E (10 MPa).

## Read first

Decisions 0019, 0032, 0035, 0037, 0038, 0050; `tests/collapse/collapse_checks.py`.

## Done when

`tests/manual/test_collapse_mui.sh` passes: run-out, top and slope checked at a = 0.5 and
2, E = 1 and 10 MPa, 10 mm.

## Log

- 2026-11-18: One session. Found and fixed material stuck against the top wall
  (ceiling counts as empty for the slab side). mu(I) + flowImpedance: 10 MPa run-out
  0.25 / 1.20 (short). With normalCap: 0.30 / 1.30 at 10 MPa, 0.40 / 1.70 at 1 MPa, all
  in band; the same rule set keeps the inclined layer passing and the tank at 73.6 /
  61.0 / 58.8 kg/s/m. Next: tutorial to the T0 rule set.
