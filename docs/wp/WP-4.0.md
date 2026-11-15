# WP 4.0: Investigate Lagrangian-Eulerian contact schemes (LS-DYNA, Abaqus, others)

State: not started. Scheduled before WP 4.1 (first package of the Lagrangian-body phase).
Branch: `wp-4.0-contact-investigation`

Requested by Sia (2026-11-15): LS-DYNA is said to handle this contact better than
Abaqus; investigate thoroughly when the time comes, before building ours.

## Goal

A recommendation, backed by sources and by small numerical experiments, for how a moving
Lagrangian body exchanges force with the Eulerian material in `celCoupling`.

Scope note: tank and silo walls are mesh boundaries with Coulomb friction
(`docs/decisions/0020`), not contact. This package concerns bodies that move through
the Eulerian mesh (footing, pile, cone, projectile: targets T1 to T3).

## Read first

- `docs/decisions/0001`, `0013`, `0014`, `0016`; `docs/reuse-survey.md` (Staubach and
  Martinelli 2023: Abaqus CEL tip stresses oscillate strongly in penetration)
- Known issues in `docs/STATUS.md` on impulsive start-up overshoot (decisions 0008, 0013)
- The void treatment decisions (0021 onward): a body displaces material much as void does

## Build

`docs/research/contact-schemes.md`, every claim tied to a source that was opened
(manual section, theory guide, paper), anything else marked "not verified":

1. **LS-DYNA.** `*CONSTRAINED_LAGRANGE_IN_SOLID`: constraint-based coupling (conserves
   momentum, not energy) against penalty-based coupling (coupling points on the
   Lagrangian segments, penalty spring on penetration relative to the tracked material);
   coupling direction options; leakage control (`ILEAK`), number of coupling points
   (`NQUAD`), penalty factor and damping, friction; the newer structured-ALE coupling
   (`*ALE_STRUCTURED_FSI`). What each is recommended for, and its known failure modes.
2. **Abaqus/Explicit CEL.** Eulerian-Lagrangian contact inside general contact: penalty
   enforcement between the Lagrangian surface and the Eulerian material interface
   reconstructed from volume fractions; requirement on mesh size ratio; documented
   leakage and penetration behaviour; contact noise in penetration problems.
3. **Other families** worth comparing: immersed boundary (direct forcing),
   cut-cell / embedded boundary with a Riemann solver at the body surface (natural for
   our cell-centred Godunov scheme: a moving wall is a face with prescribed velocity),
   fictitious domain with distributed Lagrange multipliers, and the MPM-style contact.
4. **Criteria**, applied to each: leakage, momentum and energy conservation, noise in
   the contact force (the Abaqus tip-stress problem), friction (stick and slip),
   time-step restrictions, behaviour when the body edge cuts a cell, cost, parallel
   complexity, and fit with our scheme (cell-centred, unsplit, face Riemann solver,
   `phiL` face flux, void closure).
5. **Experiments** in the 2D Python prototype (small, each under a few minutes):
   a rigid flat punch pushed into material at constant velocity, and a rigid block
   translating through material at rest, implemented with at least two candidate
   schemes (a penalty coupling and a moving-wall/cut-cell Riemann coupling). Measure
   leakage (mass inside the body), force noise against the analytic or converged
   force, and sensitivity to the penalty stiffness and to where the body edge sits
   in a cell.
6. **Whether the claim holds**: say plainly, with the evidence, whether and in what
   respect LS-DYNA's coupling is better than Abaqus's, and what that means for us.

## Done when

- `docs/research/contact-schemes.md` exists with sources (title, authors or vendor,
  year or version, link, "opened on" date), the comparison table, and the experiment
  numbers.
- `tests/test_proto_contact_schemes.sh` runs the experiments in CI and asserts the
  recommended scheme's leakage and force-noise numbers.
- `docs/decisions/NNNN-contact-scheme.md` records the recommendation and what it changes
  in WP 4.1 to 4.4 (whose rows and briefs are then updated).
- Sia is shown the recommendation before WP 4.2 and 4.3 are built on it.

## Out of scope

- The production C++ coupling (WP 4.2, 4.3). Deformable bodies (Phase 5).

## Log

(none yet)
