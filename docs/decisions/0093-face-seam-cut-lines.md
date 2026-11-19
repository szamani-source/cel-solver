# 0093: A surface seam on a mesh face keeps its cut line; the sliver fix stays parked

Date: 2026-11-19 (accuracy package A.3, thread "Fix body contact above surface")
Status: geometry fix accepted (A.3a); the sliver fix was merged in A.3b (decision 0094).

## Context

Decision 0092 parked a fix (empty slivers under a body stay empty, their full target
stays full) because steel hex8 bodies then packed the clay to 1.056 rho0 in
`test_projectile.sh` (limit 1.03) and 1.050 in `test_deformable_push.sh` (limit 1.02).

## What the packing was

Tracing the steel slug cell by cell showed the jump in two cells only, the pair either
side of x = 0.15 m, in the single step the slug's base crossed into the first full row:
their open volume fell 4 % in one step against 0.2 % elsewhere. The slug's hex8 nodes
are 9.5 mm apart from x = 0.1025, so a node line sits exactly on the mesh face
x = 0.15. `cutGeometry::faceCovered` skipped every triangle with an edge lying in the
face's plane (each gave one crossing point, not two), so the cut line along the base
was lost on that face and its covered area came from the inside test, wrong by up to
the whole face. Shifting the slug 0.1 mm off the face gave 1.017 rho0, which
confirmed it. The old mixing hid this (contact started a quarter cell early, and a
partly filled group closes void first).

## Fix (merged)

`faceCovered` counts a surface vertex on the face plane as above it (h >= 0), the same
for every triangle sharing it. A triangle edge in the plane is then the cut line of
exactly one of its two triangles (the one below), never lost and never counted twice.
`test_body_cells.sh` has a `seam` case: a box whose surface strips meet on the mesh
face x = 0.5; the cells either side must read the same fraction (old code: 2.6 % apart,
covered volume 0.6 % off; now round-off).

## The sliver fix with it (not merged)

With the geometry fix, the parked patch passes:

| Check | Main | Patch + geometry fix | Limit |
| --- | --- | --- | --- |
| Footing 2.5 mm above sand at rest, psi 0 / 30 | 0.17 m/s / blows up | 8e-17 m/s both | 1e-3 m/s |
| Projectile rigid / steel: penetration | 23.84 / 23.75 mm (-26 %) | 27.52 / 27.60 mm (-15 %) | 30 % |
| Projectile steel: density | 1.019 | 1.015 | 1.03 |
| Deformable push steel: density | 1.008 | 1.007 | 1.02 |

The full CI suite then failed `test_wedge_silo.sh`: in the 3D quarter (an STL tube)
cells at the tube wall near the free surface ran away to 1e8 m/s at t = 3.5 s and the
remap's least-squares matrix went singular. The rule never fires there (the tube's
wall faces are vertical). The crash follows round-off: writing the mixed fraction as
(gm/gV)(ga/gm) instead of ga/gV (equal in exact arithmetic) is enough to switch it on
or off. The steel push behaves the same way: with ga/gV it spikes to 1.032 rho0 at
t = 0.12 s in a corner cell beside the footing (alpha 0.9999, so on the partly filled
rules); at maxCo 0.399 or 0.401 it reads 1.014. Main passes both tests on the same
knife edge. A gap-based contact onset for slivers was tried too and was not needed
(1.026 against 1.015 without it).

## Decision

1. Merge the geometry fix and the `seam` check (A.3a). It changes no other result:
   footing, piles, STL walls, push and silo, free body, soft push, MPI, restart and
   the wedge silo pass; the steel push reads 1.016 rho0.
2. Keep the sliver fix parked. Merging it on a variant that happens to pass would
   hide the fragility. A.3b: make the free-surface small cells robust to round-off
   (wall slivers at a free surface, nearly full corner cells), then merge the patch.

## Consequences

- Any body whose surface vertices land exactly on mesh faces (structured hex8 boxes,
  STLs built on the mesh's grid) now gets exact cut fractions.
- The projectile gap stays at -26 % until A.3b; about 40 % of it is the early contact.
- `test_wedge_silo.sh` and `test_deformable_push.sh` pass on main by a margin that
  round-off can remove; a change near small cut cells may trip them for that reason.
