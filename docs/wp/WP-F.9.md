# WP F.9: Restart, MPI, tutorial, manual, v2.0.0

State: in progress (F.9a, F.9b done; F.9c restart, MPI, cost, tutorials in progress)
Branch: assigned by the thread

## Goal

Release version 2.0.0 with the gas model.

## Read first

- `docs/decisions/0125-gas-coupling-plan.md`
- `docs/spec/gas-references.md`, `tools/gas_reference.py`
- Briefs of the earlier F packages (their logs)

## Build

- Restart of p and U_g, MPI checks, `tutorials/gasHopper`, manual chapter, VERSION, release notes, CHANGELOG

## Done when

- `tests/test_gas_parallel.sh`: restart and 2/4 processors match serial to 1e-12 on the F.5 and F.8 cases; `tests/test_manual.sh` passes with the gas chapter. Speed (`tests/manual/test_gas_scaling.sh`, F.5 fine-powder hopper and the Couto hopper, 1/2/4/8 processors, 200 steps): gas on costs at most 1.5 times gas off per step, and 4 processors run at least 2.5 times faster than serial (v1 baseline: 2.6 at 620 cells per processor, decision in WP 6.3a).

Tests use the `ergunWenYu` drag law (Sia, 2026-11-24).

Validation pack: `tests/manual/test_gas_validation_pack.sh` runs the four pack tests (F.2, F.7, F.3, F.5) on 1, 2 and 4 processors, checks each against its reference and the processor counts against each other, and writes one table of results and run times to `docs/results/gas-validation.md`.

## Out of scope

- Temperature, liquids (out of v2).

## Log

Newest first. Each entry: date, what was done, what was learned, exact next step.

- 2026-11-27 (F.9c, decision 0141): restart state `<time>/uniform/poreGas`
  (`poreGas::writeState`/`readState`): the dam break parted by order 1 after a restart
  (frictionOnly's `Pprev` reset), now 0 difference; hopper already 1.8e-13, now 0.
  `tests/test_gas_restart_mpi.sh`. Dam break 1 x 2 cut through the column parts by
  round-off from the first write (zero effective stress at the start; serial with E
  x (1 + 1e-13) departs as fast; dry: flips of flowImpedance / fullBlend /
  overfillStress); parallelCheck 8e-16; outcome band from serial runs with E x (1 +- 1e-13,
  1e-12): run-out 4.88-5.73 h0, 2/4 processors 5.50/5.10. Scaling (8,704-cell hopper):
  4 processors 2.96x, gas 1.94x dry (pTolerance 1e-9 Pa, 11 solves per step); 1.50x at
  1e-4 Pa. Tutorials `tutorials/solidFluid/airImpededHopper`, `fluidisedDamBreak`,
  `tests/test_gas_tutorials.sh`. Asked Sia: default pTolerance; outcome comparison of
  the dam break. Next: apply the answers, then manual chapter, validation pack, release.

- 2026-11-24: brief written in WP F.0. Next: start when the package before it has merged.
