# 0048 STL cut-cell walls with the free surface and outlets (WP 4.3c)

Date: 2026-11-18. Status: accepted.

## Context

Decision 0047 (WP 4.3a) put fixed STL bodies into the Lagrange step and remap as
cut-cell walls, without `freeSurface`. The silo (WP 4.5) needs them with void (a free
top, material leaving through an outlet) and with the slab rules of decision 0031.

## Decision

1. **alpha per open volume.** With bodies, `alpha` is the material fraction of the
   cell's open volume `Vo`, and `rho` is mass over `Vo` (as in 0047). Covered cells hold
   no material; their ghost values (0047 item 3) now include alpha. Covered cells away
   from the material are empty ghosts.
2. **Slab sides.** A face closed by a body mirrors the cell, as a wall patch: closed
   faces are skipped in `slabSides` (neighbour fraction, smooth flag).
3. **Wall face.** Solved only where the material touches it (`touches()` along the
   face's largest axis component, as a staircase face); otherwise no force and the face
   moves with the body. Its volume flux `v.A` counts in the Lagrangian volume and, in the
   per-axis void closure, all along that largest axis (splitting it by components opened
   void in cells sliding along a tilted wall). Material that pulls away from a wall
   opens void in the remap; no mass enters (the 0047 leak counter stays at zero).
4. **Remap.** The slab rule works on the open area of each face and on the open length
   of the cell along each axis (open volume over the larger open area across that
   axis). Open volume replaces `V` in every fraction.
5. **Mixing groups.** A small cut cell and its partner take the group's Lagrangian
   volume ratio (their own volume change can be far off for a sliver), in the void
   closure and in the remap; a member's own alpha and mass are not clipped before the
   group is mixed (clipping lost material volume, so the group's density grew to 1.6
   rho0). Mixing averages alpha by material volume; an empty group stays empty.
6. **Empty is not covered.** "Covered" tests use the open fraction, not the mass: an
   empty open cell (m = 0) still goes through the void rules (with the mass test its
   material fraction from an earlier step stayed and created material volume).
7. **Outlets** are mesh patches (decision 0033) as before; STL bodies may close part of
   them (open area of each outlet face). `tutorials/stlSilo`: the T0 tank built from STL
   walls and floor, outlet = the mesh patch under the floor gap.

## Consequences

- `tests/test_stl_void.sh`: lake at rest in an STL tank 4e-17 m/s; a T0 block falling
  onto a 10 deg STL floor stops with no material through it; the WP 4.3a tilted tank
  without its lid, void top, against the lid version (see the WP 4.3c log).
- Without bodies every path is the old one bit for bit (checked on a falling block).
- **Corners of one cell with two closed sides** (an STL floor meeting an STL wall) get
  one oblique wall face (0047 item 2); under an anisotropic stress it is not in discrete
  equilibrium (a lake in an all-STL tank reached 0.3 m/s at its floor corners within 1000
  undamped steps; damped settling copes, as in the 0047 tilted tank).
  The lake test uses a mesh floor. A corner fix (one wall face per closed direction) is
  left for later.
- A cohesive block (c 5 kPa) landing in void is unstable with or without bodies
  (pre-existing; cohesion is WP 3.5).
