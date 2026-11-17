# WP 4.1: Rigid body: STL input and kinematics

State: done (2026-11-17)
Branch: `claude/project-thread-hilqye`

## Goal

Rigid bodies read from STL files move on prescribed paths inside celFoam, ready for
the contact search (WP 4.2) and coupling (WP 4.3a, 4.3b). No forces.

## Read first

- `docs/decisions/0001`, `0041`; `docs/wp/WP-4.0.md`

## Build

- `src/celLagrangian/rigidBody` (one body), `src/celLagrangian/rigidBodies` (all bodies of
  a case, output); celFoam calls `bodies.move()` each step and `bodies.write()` at write
  times.
- `tutorials/rigidBodyMotion` documents `constant/rigidBodies`.
- `tests/test_rigid_body.sh`, `tests/rigidBody/check_motion.py`.

## Done when

`tests/test_rigid_body.sh` passes: a fixed 45 degree wall, a cube pushed down by a
velocity ramp (table), a cube translating and spinning about an edge; centre, rotation,
velocity and omega within 1e-12 of the exact path, written STL within 1e-7, fixed wall
bit for bit, 2 processors write the same poses as serial.

## Out of scope

- Which cells the body cuts (WP 4.2), forces (WP 4.3a, 4.3b), bodies moved by forces.

## Log

- 2026-11-17: Done in one session. Errors 4e-16 (centre), 2e-16 (rotation), 3e-8 (STL,
  ascii output precision); test runs in about 1 s. Surfaces must be read from the
  undecomposed case's `constant/triSurface` in parallel (`globalPath`). Mesh question
  answered in decision 0041. Next: WP 4.0 (scheme), then 4.2.
