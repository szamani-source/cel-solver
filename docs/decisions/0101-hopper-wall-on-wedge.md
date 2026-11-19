# 0101: Conical hopper wall on the wedge as an STL body in the wedge's face planes

Date: 2026-11-19 (WP H.2, was WP 3.8d)

## Context

The conical hopper of Sia's validation case (decision 0100) needs a sloped wall on the
axisymmetric wedge (decision 0061). The wedge mesh keeps r-y faces, so the Cartesian
rule of decision 0031 forbids a mesh-conforming slope; the existing route for sloped
walls is an STL body acting as a cut-cell wall (decisions 0047, 0051). Nothing in the
cut-cell code assumed a Cartesian slab, so the question was how to describe a cone
on a wedge one cell thick.

## Decision

- The hopper is one closed STL body: the (r, y) polygon below and outside the cone
  (the cone extended below the floor and out past R) extruded across the wedge. A
  vertex (r, y) becomes (r cos t, y, +-k r sin t), t the wedge half angle. The straight
  line between the two lies at x = r cos t, the plane of the wedge cell's own flat
  faces, so the cut is the mesh's conical face exactly, for any k.
- k = 6 puts the two caps outside the mesh's hoop extent, so no body surface lies on
  a wedge face (decision 0093 had trouble with surfaces on mesh faces), and the
  polygon edges are split into 10 mm pieces so no triangle's box spans the body.
  The caps are fanned from the centroid (a fan from a corner made degenerate
  triangles). Caps must face outward: with them reversed, no cell was covered and
  the fill was a full cylinder (10 % too much mass).
- Generator: `tests/hopper/hopper_case.py` (`hopper_stl`, `couto_case`).
- Done-when, reworded: the cylinder pressure is checked against Janssen with
  K = nu/(1 - nu), the ratio an elastic fill reaches under gravity (the run's interior
  ratio is 0.425 against 0.429), not K_A = 0.338: a fill at rest is not in the active
  state P1's Nedderman formula assumes.

## Consequences

- `tests/test_hopper_wedge.sh` (about 4 minutes, 25 mm cells): mass 1501.8 kg against
  the geometry's 1501.8 kg; closed outlet holds (mass to 9e-15); cylinder pressure at
  GA1 to GA7 within -4.7 to +6.0 % of Janssen (K 0.43); transition vertical stress
  12.62 kPa, inside the experiment's 11.65 to 13.35 kPa.
- The model's cylinder pressures also sit inside or near the experiment's static
  ranges (GA1 5.3 against 4.0 to 4.5 kPa is the furthest), which K_A Janssen did not
  reach above GA3 (decision 0100).
- At rest the hopper wall reads 11.4 to 11.8 kPa, nearly uniform, against GA0's 8.6 to
  9.6 kPa. That is H.3's question.
- The same construction serves any axisymmetric wall shape (insert, cone-in-cone).
