# 0016 Re-plan around 2D silo discharge (target T0)

Date: 2026-11-14. Status: accepted (Sia named the problem on 2026-11-14).

## Context

Sia's first specific target is 2D discharge of a cohesionless bulk solid from a
rectangular tank through a bottom opening (`docs/spec/target-problems.md`, T0). The
original plan went: several materials and void (Phase 3), rigid bodies and contact
(Phase 4), deformable bodies (Phase 5).

## Decision

- T0 needs no Lagrangian body. Phases 4 and 5 (the "Lagrangian" half of CEL) are
  deferred until T0 runs and is validated. They stay in the plan unchanged.
- Phase 3 is rewritten as "silo discharge": gravity, frictional walls, one material
  plus void (free surface), an outlet that opens, and validation against Janssen and
  Beverloo. Several materials in one cell (old 3.1, 3.3, 3.4) move to a later phase 3B.
- Each new physical ingredient is tried first in the 2D Python prototype (packages
  ending in "p"), as in Phase 1, because that found the scheme problems cheaply. The
  prototype packages do not depend on the C++ Lagrange step and remap (WP 2.3, 2.4), so
  they start now.
- The void treatment (how a partly filled cell behaves) is decided in WP 3.2p with
  numbers, not here.

## Consequences

- The first end-to-end result Sia sees is a prototype silo discharge (WP 3.3p), then
  the same in `celFoam` (WP 3.4).
- The code is, for now, an Eulerian large-deformation solid solver; it becomes "coupled
  Eulerian-Lagrangian" when Phase 4 resumes.
