# 0092: Empty slivers cut off by a body: the cause is found, the fix is not kept yet

Date: 2026-11-19 (accuracy package A.2, thread "Fix accuracy issues before packaging")
Status: finding recorded; code change parked in `docs/patches/A2-empty-sliver-mixing.patch`.
The steel packing is explained and fixed in decision 0093; the patch stays parked (A.3b).

## Context

Decision 0091 found that a rigid footing held 2.5 mm above sand at rest disturbs the
sand (`tests/manual/test_body_above_surface.sh` fails), and that associated flow blows up
there before contact.

## Cause

Moving the footing's base shows where the problem is (10 mm cells, held 0.2 s):

| Base above the surface | First void row at the base | Material at rest? |
| --- | --- | --- |
| 2.5 mm | 25 % open: a small cut cell, mixed with the full cell below | no (psi 0: 0.17 m/s; psi 30: blows up) |
| 7.5 mm | 75 % open: not small, not mixed | yes (1e-17 m/s), psi 0 and 30 |
| 12.5 mm | second row | yes |

`mixSmallCells` (decision 0047) spreads each group's material evenly over the group.
When the small cell is an empty sliver, the full cell below it becomes 80 % full and
takes the partly filled free-surface rules, though its material still reaches its top
face. At the footing's corners this sets the top cells moving, and with psi = phi the
corner cell's stress grows without bound.

## What was tried

The rule was: when a group's small cells hold no material and its target is full and
not cut by the body, keep all the material in the target and leave the slivers empty
(conservative to round-off). The patch is in `docs/patches/`. Results:

- `tests/manual/test_body_above_surface.sh` passes: 8e-17 m/s for psi 0 and psi 30.
- Rigid projectile into clay: 27.5 mm against 23.75 mm before, so the gap to the
  Prandtl energy estimate (32.3 mm) falls from -26 % to -15 % (decision 0074). The old
  mixing made contact start up to a quarter cell early.
- Sand footing psi 0: no change (N_gamma 11.2 to 10.2 at 10 mm). psi 30 now reaches
  contact (q about 19 kPa, N_gamma about 21) and blows up later in the push, at 0.6 s.
- STL walls, body cells, clay footing, pile, rough pile and STL silo tests pass.
- **Regression:** the steel (hex8) bodies pack the material at impact, 1.056 rho0 in
  `test_projectile.sh` (limit 1.03; 1.13 with maxCo 0.38) and 1.050 in
  `test_deformable_push.sh` (limit 1.02). Main gives 1.019 and 1.008. So the change
  is not merged.

## Decision

1. Keep main unchanged; record the cause and the parked patch.
2. The fix needs a contact onset for slivers that does not depend on the mixing:
   the target's top face should see the body once the sliver closes, without a jump
   when a full cell is first cut. This is its own package (A.3).
3. The sand footing excess is not explained by the hovering effect (psi 0 unchanged).
   What remains: the associated blow-up during the push, and resolution (B is 10 cells).

## Consequences

- About 40 % of the projectile shortfall (known issue 3) is the early contact from
  the mixing, not the clay model.
- `tests/manual/test_body_above_surface.sh` stays failing until A.3.
