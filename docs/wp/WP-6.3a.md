# WP 6.3a: Parallel scaling

State: done (2026-11-19). Split from WP 6.3 (the report is 6.3b).
Branch: the thread's `claude/project-thread-*` branch

## Goal

Know what parallel runs buy for the T0 tank, to plan larger runs.

## Read first

- `docs/results/scaling.md`; `tests/scaling/scaling_run.py`

## Done when

`tests/manual/test_scaling.sh` runs (about 1.5 min) and its table is recorded in
`docs/results/scaling.md`: 4 processors give 2.6 (10 mm) and 3.3 (5 mm) times serial.

## Out of scope

- More than 4 cores (the cloud container has 4); 3D; Sia's laptop (final items).

## Log

- 2026-11-19: Measured. 4 us per cell per step serial; 65 to 82 % efficiency on 4 cores
  depending on cells per processor.
