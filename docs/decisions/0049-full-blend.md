# 0049 Blending the full and partly filled rules (WP 3.6h)

Date: 2026-11-18. Status: accepted (option kept, off by default; done-when not met).

## Context

Decision 0046: the T0 tank rate moves by up to 30 % with the line between "full" and
"partly filled" cells (nearly-full snap of decisions 0044 and 0045). Dilatancy did not
help (decision 0048). WP 3.6h tries the numerical alternative proposed there: replace
the hard switch with a gradual blend. Done-when: tank bands at 1 and 10 MPa within 10 %
of 0.1 MPa at 10 mm, `tests/manual/test_inclined_layer.sh` still passing.

## What was built

`fullBlend yes` (with `freeSurface yes`) and optional `fullBlendVoid (v1 v2)`, default
(1e-4 1e-2). No nearly-full snap; alpha keeps its value. A cell's weight w falls from 1
at void v1 to 0 at void v2, linear in log(void):

- face states: w times the second-order state plus (1 - w) times the first-order one
  (w of a cell = smallest of its own and its neighbours' weights);
- Hancock stress predictor times w;
- flow viscosity (flow impedance) times the cell's own w;
- slab geometry: cells with w >= 1/2 (void below sqrt(v1 v2)) stand upright and touch
  all faces; the remap keeps its 0.5 threshold on the weight.

Without `fullBlend` the weights are 0 or 1 and every result is unchanged bit for bit
(the code takes the same branches).

## Results

Tank bands, 10 mm, three runs each, kg/s per m (`tests/manual/test_discharge_scatter.sh`
with `CEL_EXTRA_MATERIAL='fullBlend yes; ...'`); inclined layer profile error, 5 mm:

| Variant | 0.1 MPa | 1 MPa | 10 MPa | Layer 0.1 / 1 / 10 MPa |
| --- | --- | --- | --- | --- |
| Snap 1e-4 with density bound (decisions 0045, 0046) | 65.8 | 52.0 | 60.8 | 3.6 / 0.9 / 0.8 % pass |
| No snap (decision 0046) | 65.5 | 60.5 | 63.1 | fails |
| Blend 1e-4 to 1e-2 | 68.7 (68.6 to 68.9) | 53.7 (52.3 to 56.0) | 41.7 (41.3 to 42.2) | 3.6 / 0.9 / 0.9 % pass |
| Blend 1e-6 to 1e-4 | 62.5 (61.8 to 63.2) | 42.6 (42.5 to 42.7) | 47.0 (46.5 to 47.6) | 3.6 / 1.3 / 944 % fail |

## Findings

1. The blend cuts the run-to-run scatter (bands 0 to 7 %, 10 MPa 2 % instead of 18 %),
   so the switch's discontinuity did cause the scatter.
2. It does not remove the E dependence. The wide blend keeps the layer but slows
   10 MPa by a third; the narrow one slows 1 MPa and lets the stiff layer's top cell
   slide off. Together with decision 0046's table, six treatments of loosened cells
   give 1 MPa rates from 42 to 67 kg/s/m with no trend that a setting could tune.
3. So the E dependence is in how the partly filled rules themselves treat a loosening
   stream (first-order states, no stress predictor, stress-free opening), not in where
   or how sharply they switch on. Neither a material law (0048) nor the switch (this
   note) fixes it.

## Decision

- Keep `fullBlend` as an option, off by default (snap stays the default: it is the
  only setting that passes the layer at all three E with a 10 MPa rate within 10 %).
- Close WP 3.6h with its done-when not met. Stop tuning the void rules for E
  independence. The 1 MPa dip (about 20 %) stays a known issue.
- Proposed next, for the coordinator and Sia: measure the tank at 5 mm and at T0's
  real stiffness only, and decide whether a rate uncertainty of about 20 % across
  artificial stiffness is acceptable for now, before more work on partly filled cells.

## Consequences

- New optional entries `fullBlend`, `fullBlendVoid` (documented in `eulerianSolid.H` and
  here; no tutorial uses them).
