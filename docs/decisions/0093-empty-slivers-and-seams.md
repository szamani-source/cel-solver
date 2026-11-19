# 0093: Empty slivers stay empty; a surface seam on a mesh face keeps its cut line

Date: 2026-11-19 (accuracy package A.3, thread "Fix body contact above surface")
Status: accepted; replaces the parked patch of decision 0092.

## Context

Decision 0092 found that small-cell mixing spreads material into the empty sliver
under a body just above the surface, and parked a fix (keep empty slivers empty and
their full target full) because steel hex8 bodies then packed the clay to 1.056 rho0
in `test_projectile.sh` (limit 1.03) and 1.050 in `test_deformable_push.sh` (limit
1.02). The guess was a jump in contact when a full cell is first cut.

## What the packing was

Tracing the steel slug cell by cell showed the jump in two cells only, the pair either
side of x = 0.15 m, in the single step the slug's base crossed into the first full row:
their open volume fell 4 % in one step against 0.2 % elsewhere. The slug's hex8 nodes
are 9.5 mm apart from x = 0.1025, so a node line sits exactly on the mesh face
x = 0.15. `cutGeometry::faceCovered` skipped every triangle with an edge lying in the
face's plane (each such triangle gave one crossing point, not two), so the cut line
along the base was lost on that face and its covered area came from the inside test:
wrong by up to the whole face. The old mixing hid this (contact started a quarter cell
early and the partly filled group closed void first). Shifting the slug 0.1 mm off the
face gave 1.017 rho0 with the parked patch, which confirmed it.

A gap-based contact onset for slivers (the sliver's wall face acting on its full target
once the body closes the gap within the step) was tried too. It passes the tests but
is not needed: without it the packing is 1.015 rho0 (steel) and 1.013 (rigid); with it
1.026, from a sharper impact wave. Not kept.

## Decision

1. Merge the empty-sliver rule of decision 0092 (`mixSmallCells`): when a group's
   small cells hold no material and its target is full and not cut by a body, the
   slivers stay empty and the target keeps all the mass.
2. `faceCovered` counts a surface vertex on the face plane as above it (h >= 0), the
   same for every triangle that shares it. A triangle edge in the plane is then the cut
   line of exactly one of its two triangles (the one below), so it is never lost and
   never counted twice. Exact to round-off, as before.
3. `tests/test_body_above_surface.sh` moves from `tests/manual` to CI. The body-cell
   test gets a `seam` case: a box whose surface strips meet on a mesh face; the cells
   either side read the same fraction (the old code: 2.6 % apart, volume 0.6 % off).

## Results (10 mm cells)

| Check | Main before | Now | Limit |
| --- | --- | --- | --- |
| Footing 2.5 mm above sand at rest, psi 0 / 30 | 0.17 m/s / blows up | 8e-17 m/s both | 1e-3 m/s |
| Projectile, rigid: penetration, density | 23.84 mm, - | 27.52 mm (-14.9 %), 1.013 | 30 %, 1.03 |
| Projectile, steel hex8: penetration, density | 23.75 mm, 1.019 | 27.60 mm (-14.6 %), 1.015 | 30 %, 1.03 |
| Steel against rigid | -0.3 % | +0.3 % | 5 % |
| Deformable push, density steel / rigid | 1.008 | 1.007 / 1.009 | 1.02 |

## Consequences

- The projectile gap to the Prandtl energy estimate falls from -26 % to -15 %
  (known issue 3): contact no longer starts a quarter cell early.
- Contact with a body now starts when it reaches the material, so the impact is
  sharper; the clay's peak density comes from the elastic wave (steel 1.015 rho0).
- The sand footing excess (known issue 2) is not changed at psi 0 (decision 0092);
  the next package looks at psi 30 during the push.
