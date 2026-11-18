# 0042: Couple bodies and STL walls as cut-cell moving walls, not penalty springs

Date: 2026-11-17 (WP 4.0). Shown to Sia before WP 4.3a and 4.3b are built.

## Context

Sia asked (2026-11-15) for a thorough look at how LS-DYNA, Abaqus and others couple
Lagrangian bodies to Eulerian material before ours is built, and (2026-11-17) for STL
walls (decision 0041). Sources and experiments: `docs/research/contact-schemes.md`.

LS-DYNA (CLIS, S-ALE) and Abaqus CEL both use penalty coupling. Both document leakage
and offer remedies. In the 2D prototype (`prototype/contact2d.py`) a penalty coupling
showed the usual trade: soft springs leak (1 %) and under-predict the force by 15 %;
stiff springs ring (contact stress noise 27 % at k = M/h, over 100 % at 10 M/h) and get
the friction of a sliding block wrong by 11 to 95 %. A cut-cell coupling, where the body
surface is a face with the body's velocity and our acoustic Riemann solver gives the
traction, gave the exact contact stress within 0.6 % with 0.2 % noise whatever the
plate's position in the cell, no material through the wall, and the sliding
acceleration on a 30 deg frictional wall within 0.3 %.

## Decision

1. Bodies and STL walls act on the material through cut cells: faces are cut to their
   open area, and the body surface inside a cut cell is one more face whose state is the
   one-sided Riemann solution with the body's velocity (no tension: the face separates;
   Coulomb cap on the shear with the wall friction). No mass crosses the body. The force
   on the body is minus the sum of these face forces.
2. Small cut cells (open fraction below 0.5) are mixed with their neighbour away from the
   body (velocity and stress, mass and momentum kept); a cell the body covers hands its
   content to that neighbour. A donor never gives more mass than it holds.
3. No penalty coupling is built. If deformable bodies (Phase 5) need it, revisit.

## Consequences

- WP 4.2 must give, per cell and step: body fraction, open fraction of each face, the
  wall segment (area times normal, taken as minus the sum of open face areas times their
  normals so a uniform stress gives no force), the body velocity there, and the mixing
  neighbour. Its brief is updated.
- WP 4.3a puts the wall face into the C++ Lagrange step next to the box walls (decision
  0030) and the remap (no flux through covered face parts); the hydrostatic gravity
  reconstruction and the void rules must see the cut faces.
- The prototype is first order, linear elastic and has rectangles and half-planes only;
  the C++ version meets curved STL facets, second-order face states and the plastic
  models. The mixing step is the part most likely to need work (noise 6 % on the block
  with its edge on a face, 2 % off it).
- Thin walls thinner than a cell need a closed STL with thickness of at least a cell, as
  Abaqus advises for its own coupling (S5); recorded for WP 4.2.

## Addendum (2026-11-17): Sia's note on LS-DYNA

Sia: LS-DYNA handles contact points better than Abaqus; use its method or a similar
one. We built an LS-DYNA-style coupling-point variant in the prototype
(`docs/research/contact-schemes.md`, follow-up; manual details not verified). It is
smoother than a cell-overlap spring (Abaqus-like), which supports Sia's point, but
it does not beat the cut-cell wall on any measure; with leak control it becomes as
noisy as the plain spring. The decision stands: cut-cell moving walls. What we take
from LS-DYNA into them:

1. **Coupling-point density**: the wall face in each cut cell comes from exact clipping
   of the STL facets, the limit of many coupling points per cell. WP 4.2 warns when an
   STL feature is thinner than about two cells, and when facets are much larger than a
   cell on a curved surface (S3: mesh ratio near 1:1, a few cells around thin parts).
2. **Leakage control**: per body, report the material found inside it (covered cells,
   overfilled cut cells) at every write, as LS-DYNA's leak diagnostics; such material is
   handed out of the body (point 2 above). If WP 4.3a still finds leakage, a soft
   penalty on it, like ILEAK, is the fallback, off by default.
3. **Damping and stiffness**: nothing to tune. The wall Riemann solution already acts
   with the material's impedance, the exact non-reflecting damping; no PFAC.
4. **Compression only and friction** (as DIREC and the friction coefficient): no tension
   across the wall, Coulomb cap on the shear, per body friction in `constant/rigidBodies`.
5. **Keep the comparison runnable**: the `points` and `penalty` schemes stay in
   `prototype/contact2d.py`, so a penalty coupling can be revisited for deformable bodies
   (Phase 5), where LS-DYNA's approach is the usual one.
