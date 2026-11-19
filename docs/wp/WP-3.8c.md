# WP 3.8c: Cylindrical silo on the wedge

State: done
Branch: `claude/project-thread-53vkqx`

## Goal

A cylindrical silo with a central circular outlet discharging on an axisymmetric
wedge, instead of a 3D quarter model (Sia, 2026-11-18).

## Read first

- `docs/decisions/0061-axisymmetric-wedge.md`, `docs/decisions/0062-wedge-silo-outlet-resolution.md`

## Build

- No solver changes: the outlet patch sits on the axis block of the wedge mesh.
- `tests/wedgesilo/wedge_silo_checks.py`: wedge mesh with two radial blocks; quarter
  mesh with `symX`/`symZ` planes and an STL tube wall (`rigidBody`, 20 deg).

## Done when

- `tests/test_wedge_silo.sh` passes: closed outlet holds, mass balance 1e-10; 3D
  Beverloo C (k 1.5, d 1 mm) within 0.45 to 0.75 for both models; wedge C within
  10 % of the quarter model's.

## Out of scope

- STL walls on a wedge (3.8d); off-centre outlets (not axisymmetric).
