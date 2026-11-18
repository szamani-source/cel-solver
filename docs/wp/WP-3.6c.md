# WP 3.6c: flow impedance (less numerical resistance in flowing material)

State: closed, done-when not met at 5 mm (2026-11-17); rest split to WP 3.6c2
Branch: `claude/project-thread-whja6i`

## Goal

Reduce the numerical dissipation at faces where the material flows, and make the shear
cap use the mu(I) friction, so the inclined layer of WP 3.6b follows the Bagnold profile
at every E.

## Read first

- Decisions 0036, 0039, 0040; `docs/wp/WP-3.6b.md`
- `tests/layer/inclined_layer.py`; `src/celEulerian/eulerianSolid/eulerianSolidLagrange.C`
  (`flowViscosity`, `flowJump`, the Hancock predictor, the face loops)

## Done when

The inclined layer within 5 % of the exact solution (surface speed and profile) at
E = 0.1, 1 and 10 MPa, 5 mm; `tests/test_collapse.sh`, `tests/test_outlet.sh`,
`tests/test_discharge.sh` still pass. Test: `tests/manual/test_inclined_layer.sh`
(about 6 minutes; fails at 5 mm). CI part: `tests/test_inclined_layer.sh` (guard).

## Out of scope

- Shear-band mesh dependence (decision 0038), the silo with mu(I) (WP 3.6d).

## Log

- 2026-11-17: One session. Scaled each face term on the layer: three sources of
  resistance (decision 0040): the elastic Hancock predictor (largest with the cap on),
  the elastic shear jump (largest without the cap), the static-phi cap (lets faces slip,
  the layer then runs away). Built `flowImpedance yes`: mu(I) tangent viscosity in the
  face shear jump between flowing cells and at walls, plastic predictor along the flow
  direction, cap skipped between flowing cells and at phi_m elsewhere, 3 % band for
  "flowing". Body speeds now 0.65 to 0.97 of exact at 0.1 MPa (was 0.14 to 0.20), 0.41 to
  0.92 at 1 MPa, 0.03 to 0.74 at 10 MPa; 2.5 mm passes 5 % at 0.1 and 1 MPa. Not met:
  5 mm (7.8 % at 0.1 MPa, floor cells), the 10 MPa body (floor cell nearly stuck), the
  top cell sliding at 1 and 10 MPa (stress-free at 10 MPa). Silo, collapse and outlet
  tests pass (switch off there). Next: WP 3.6c2 (10 MPa floor and top cell), then 3.6d.
