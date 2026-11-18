# WP 3.6h: blend of the full and partly filled rules

State: closed, done-when not met (2026-11-18)
Branch: `claude/project-thread-4zdz5u` (with WP 3.6g, PR #44)

## Goal

Replace the nearly-full snap with a gradual blend so the T0 tank rate agrees across E
and the inclined layer still passes (decision 0048 proposal).

## Read first

Decisions 0044, 0045, 0046, 0049; `slabSides` and `fullWeight` in
`eulerianSolidLagrange.C`.

## Done when

`tests/manual/test_discharge_scatter.sh` with `fullBlend yes`: 1 and 10 MPa bands within
10 % of 0.1 MPa at 10 mm; `tests/manual/test_inclined_layer.sh` passes.

## Log

- 2026-11-18: One session. Built `fullBlend` / `fullBlendVoid`. Blend 1e-4 to 1e-2:
  68.7 / 53.7 / 41.7 kg/s/m, layer passes. Blend 1e-6 to 1e-4: 62.5 / 42.6 / 47.0, layer
  fails at 10 MPa. Scatter falls to 0 to 7 %, E dependence stays. Stopped tuning;
  next is a decision with Sia (decision 0049).
