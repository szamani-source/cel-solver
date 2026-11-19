# WP 6.2a: Restart

State: done (2026-11-18). Split from WP 6.2 (mass scaling to 6.2b).
Branch: the thread's `claude/project-thread-*` branch

## Goal

A celFoam run restarted from any time it wrote continues exactly as the run that wrote
it, with bodies, outlets and inlets, serial and parallel.

## Read first

- `docs/decisions/0081`; `applications/solvers/celFoam/celFoam.C` (time step);
  `rigidBody::writeState`, `hex8Body::writeState`, `eulerianSolid::writeState`

## Build

- `<time>/uniform/eulerianSolid` and `<time>/uniform/rigidBodies` written at each
  write time (every processor), read at the start time.
- celFoam's own step: stable, landing on write times counted from zero and on endTime.
- postProcessing files appended on a restart.

## Done when

`tests/test_restart.sh` passes (about 40 s): T0 tank restarted while discharging; tank
with a free block and a hex8 footing in contact, serial and on 2 processors; fields,
state and bodyForces.dat after the restart identical bit for bit.

## How to restart a case

Keep `writeFormat binary` (or `writePrecision 17`) for exact restarts, then set
`startFrom latestTime;` (or `startTime <t>`) in system/controlDict and run celFoam
again, in parallel from the processor directories as usual.

## Out of scope

- Mass scaling (WP 6.2b).
- Restart of celVoidCheck/celSiloCheck test drivers (they start from their fill).

## Log

- 2026-11-18: With OpenFOAM's adjustable step the restarted run took steps differing in
  the last digit (write time counted from the start time) and parted. celFoam now
  computes the step (decision 0081). Also found: free bodies refused to start after
  t = 0; now they need their state file instead. All three cases identical bit for bit.
