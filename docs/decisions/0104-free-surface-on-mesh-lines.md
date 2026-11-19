# 0104: Free surface along mesh lines on angled cells

Date: 2026-11-19 (WP G.2)

## Context

The free-surface rules (decision 0031: slab sides, touching faces, void
closure, swept-slice remap, layer spreading) work per x, y, z axis and stopped
with an error on any face off the axes. Sia's wall-aligned hopper meshes
(decision 0103) have leaning cells.

## Decision

1. Each cell of an angled mesh gets three mesh lines. Its faces pair up by most
   nearly opposite normals; a pair's line joins the two face centres; a face
   left alone (the outer face of a wedge axis prism) gives the line from the
   cell centre. For a parallelepiped a face's area vector then has no component
   along the other two lines, so it is classified exactly.
2. The lines are labelled by the axis each lies closest to (the permutation with
   the largest sum of |t_a . e_a|) and point along +axis, so neighbouring cells
   of a structured mesh agree on labels. Meshes whose lines turn more than about
   45 degrees from the axes (O-grids) are not covered.
3. Every rule that took a face's axis from its area vector S now takes it from
   S in the cell's line frame (`inLines`): slab sides and ceilings, touching
   faces, the void closure's volume per axis, the remap's swept slice, and the
   rule 7 weight along gravity. Cell lengths h_a = V/A_a with A_a half the
   area of the faces on line a (perpendicular height of the pair).
4. A coupled (processor) neighbour's face is classified in the local cell's
   frame: the two agree on a structured mesh.
5. Rule 9 (layer spreading) compared alpha differences along the mesh lines.
   On leaning columns the difference "along gravity" picked up the surface
   slope, and cells leaning against the flow spread a thin tongue along the
   floor (run-out 0.40 m against 0.26 on the Cartesian mesh). On angled
   meshes the differences are turned into the gradient's x, y, z components
   first (D_a ~ 2 h_a t_a.grad alpha); the tongue goes (0.28 m).
6. Axis-aligned meshes take none of these paths (no frame), so their results are
   unchanged bit for bit. STL walls with a free surface on an angled mesh stop
   with an error until WP G.5.

## Consequences

- `tests/test_angled_free_surface.sh` (CI, about 75 s): lake at rest on leaning
  cells (mass exact, flat surface, s_yy within 4.5 % at the block kink, 0.17 %
  on average); column collapse on cells leaning both ways against Cartesian:
  front within 2.7 cells, height equal, profile RMS 0.7 to 0.9 cells.
- The largest lake stress error sits where the mesh lines kink between the two
  blocks; smoother meshes (graded blocks) should do better.
- Walls, outlet, inlet on angled faces (G.3) and the wedge axis prisms (G.4)
  come next.
- Tried: the void closure's normal stress and strain along the unit mesh line
  (t.s.t, adding the strain difference times t t) instead of the x, y, z
  components. The leaning-left collapse grew its tongue again (front 0.44 m)
  and the lake's largest stress error rose to 6.3 %, so the closure stays on
  the x, y, z components.
