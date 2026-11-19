# WP 6.1: MPI with Lagrangian bodies

State: done (2026-11-18). Needs WP 4.3d, 4.5, 5.3a.
Branch: the thread's `claude/project-thread-*` branch

## Goal

Runs with rigid, free and deformable bodies and STL walls give the same results on N
processors as serial.

## Read first

- `docs/decisions/0080` (what "agree" means), `0046` (bands of runs), `0051` (mixing)
- `tests/mpi/mpi_checks.py`

## Build

- `tests/test_mpi_bodies.sh` and `tests/mpi/mpi_checks.py`. No solver change was needed.

## Done when

`tests/test_mpi_bodies.sh` passes (about 8 min): free block drop and hex8 footing push
on 2 processors (split through the body), STL tank with an outlet on 4 (2 x 2); first
100 steps match serial to 1e-12 (or 10 times a perturbed serial run's difference);
force at rest within 2 %, footing force within 5 %, discharge rate within 5 %; mass,
hand-over loss and leak as serial.

## Out of scope

- Small cut cells mixing across processor boundaries (not needed so far; reported as
  unmixed if it happens).
- Speed-up measurements (WP 6.3a).

## Log

- 2026-11-18: First runs part from serial after 222 steps (push) or 494 (drop); a serial
  run with a 4e-14 perturbation parts the same way at the same step, so it is round-off
  through the switches, not a parallel error (decision 0080). Results: drop force /
  weight 1.0000 both (perturbed 1.0002), height 0.006 mm apart; footing 37.8 against
  39.3 N (perturbed 38.7); tank rate 68.26 serial, 4 processors and perturbed. A
  stress-free fill differs from step 2 in the cap counts with any perturbation, hence
  the relative exact-phase bound. Lesson: after changing a class layout, rebuild the
  test drivers (celVoidCheck ran with an old eulerianSolid layout and smashed its stack).
