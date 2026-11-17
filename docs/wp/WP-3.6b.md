# WP 3.6b: steady inclined layer (Bagnold profile) with mu(I)

State: closed, done-when not met (2026-11-17); measurement recorded (decision 0039)
Branch: `claude/project-thread-747gra`

## Goal

Measure the numerical viscosity of celFoam's face Riemann solver on plastic flow with a
steady inclined layer, mu(I) friction (decision 0037), E = 0.1, 1, 10 MPa.

## Read first

- Decisions 0036, 0037, 0039; `docs/rate-dependent-friction-survey.md`
- `tests/layer/inclined_layer.py` (set-up and exact solution in its docstring)

## Done when

Surface speed and profile within 5 % of the exact solution at 0.1 MPa; the error
against E measures the numerical viscosity (recorded either way).
Test: `tests/manual/test_inclined_layer.sh` (about 3 min; fails now). CI part:
`tests/test_inclined_layer.sh` (about 10 s).

## Out of scope

- The plastic impedance at flowing faces (WP 3.6c), shear-band regularisation (decision 0038).

## Log

- 2026-11-17: One session. Layer 0.05 m on 31 deg, x periodic, gravity tilted on an
  untilted mesh, sticking floor, void above or a smooth lid, from the exact steady state.
  First runs slid uphill: the hydrostatic face reconstruction assumed a column along g;
  fixed for tilted g only (decision 0039 item 1; at rest at 25 deg: 2.4e-5 m/s, was
  4.8e-4). Result: the body flows at 0.14 to 0.20 of the exact speed at 0.1 MPa, 0.05 to
  0.09 at 1 MPa, 0.01 to 0.02 at 10 MPa (5 mm); 2.5 mm doubles it at 0.1 MPa. Numerical
  viscosity about 0.06 rho c_s h (2.8 / 6.5 / 29 Pa s), against a mu(I) tangent
  viscosity of 0.8 Pa s. The top cell slides over the body at 1 and 10 MPa (likely the shear
  cap at static phi; not checked). Next: WP 3.6c.
