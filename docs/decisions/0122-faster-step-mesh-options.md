# 0122: Optional mesh options for a longer time step (off by default)

Date: 2026-11-23. Thread "Build case generator and mesh checker" (WP G.6).

## Context

celCheck (S.3) showed what sets the step on the Couto meshes: on the standard mesh the
axis cell in the spout (its hoop length is its distance from the axis), on the fine
mesh the 1.21 mm columns at the outlet edge, squeezed by the cylinder's wall grading
that the hopper block shares. Sia (2026-11-23): accuracy first; speed is a possible
later step. So the options exist but are off.

## Decision

- conicalHopper `mesh: outletColumns: even` makes the hopper block's columns uniform
  across the outlet (edgeGrading: the wall grading stays on the shell side), and
  `axisColumn: f` makes the column on the axis f outlet columns wide. flatBottom has
  `axisColumn` too. Both advanced, defaults shared / 1 (the validated mesh, unchanged
  bit for bit). blockMesh only.
- Gains measured by celCheck, equal to celFoam's step:

| Mesh | Default | even, axis 1 | even, axis 2 | even, axis 3 |
| --- | --- | --- | --- | --- |
| Couto standard (7,584) | 5.85e-6 s | 4.59e-6 | 9.03e-6 (1.54x) | 1.01e-5 (1.72x) |
| Couto fine (20,212) | 3.75e-6 s | 3.50e-6 | 6.93e-6 (1.85x) | 7.54e-6 (2.01x) |
| flatBottom silo example | 1.11e-5 s | | 1.88e-5 (1.69x) | |

  "even" alone is shorter (the axis column gets thinner); it needs the wider axis
  column. With axis 3 the outlet-edge cells set the step again.

## Consequences

- Not validated: whether the wider axis column changes the at-rest pressures or the
  discharge needs a settling and discharge run on the Couto case. Until Sia asks for
  that run, these stay off and the manual says so.
- `tests/test_celsetup_fastdt.sh`: defaults unchanged, the gains, Mesh OK, mass placed,
  celFoam runs at the predicted step.
