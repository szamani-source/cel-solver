# 0106: Free surface on slopes meshed with angled cells (WP G.3, in progress)

## Context

G.1 and G.2 tested angled cells where one family of mesh lines stays
horizontal (leaning columns, the wall-aligned hopper). G.3 asks for walls
and free surfaces on faces at any angle. Its test (`tests/manual/test_angled_walls.sh` (moves to `tests/` when it passes))
slides an elastic layer 60 mm thick down a 30 degree frictional slope
(periodic strip, void above) on three meshes: Cartesian with gravity
turned (reference), the strip turned 30 degrees (square cells, no face on an
axis) and a sheared strip (lines along the slope and vertical lines).

## Findings (2026-11-20)

- Before any change, the turned strip blew up within 0.1 s (87 % of the mass
  lost through the top). The void closure (decision 0021 rule 5) took the
  normal stress and strain of each line from the x, y, z components.
- Decision: in cells whose mesh lines are square to each other but turned,
  the void closure works in the lines' frame (`lineRot_`, the stress and
  strain turned into it and back). Cells with leaning lines keep x, y, z,
  as decision 0104 found that better there; axis-aligned meshes take no new
  path. G.1, G.2 and Cartesian results are unchanged.
- Result: the mean speed down the slope is now within 0.1 % (turned) and
  0.7 % (sheared) of the closed form, and the Cartesian reference within
  0.01 %. Still failing: the layer drifts off the wall (normal speed 4 to
  7 % of the slide), a little mass leaves through the top (1e-4 turned,
  3e-6 sheared), and a layer that should rest at friction 60 deg creeps at
  0.3 to 0.5 m/s. On the sheared strip the row against the wall runs faster
  than the rows above it.
- Tried and dropped: the hydrostatic jump of a layer normal to the two mesh
  lines across gravity (instead of the x, y, z layer of decision 0039).
  Correct in principle for a slope, but it made both meshes worse (turned
  slide +13 %, sheared +1.8 % with more lift-off).
- Without a free surface (celFoam, full strip under a frictionless lid) the
  turned and sheared strips slide at the closed-form speed in every cell
  (0.07 %): walls and friction on angled faces are right; what remains is
  in the free-surface rules. The test checks this first.
- Also tried and dropped: limiting the face states in the turned cells'
  line frame together with the layer jump above. The cell at rest improved
  (creep 0.15 against 0.52 m/s) but the slide got worse (+14 %); the turned
  lake still drifts from the Cartesian one within the first steps (shear
  appears in the top full row), so another rule is still frame-dependent.
- A full slab under a frictionless lid sticks even on the Cartesian mesh
  (wall row sheared against the lid): a separate free-surface-under-a-lid
  issue, noted for later.

## Consequences

- The hopper and Couto cases (horizontal lines, gravity along a line family)
  are not affected.
- G.3 stays open: next suspects are the wall-row face states (mirror ghost
  plus the x, y, z hydrostatic jump on faces not normal to gravity) and the
  surface rules for faces at an angle to gravity.
