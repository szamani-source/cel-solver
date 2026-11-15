# WP 3.3r: Prototype discharge rate at T0 stiffness (split from WP 3.3q)

State: not started (split 2026-11-15, decision 0027)
Branch: thread branch

## Goal

The T0 discharge rate in `prototype/discharge2d.py` does not depend on E between 0.1
and 10 MPa (T0's own value), so the prototype can be run, and later the C++ code
checked, at the real stiffness.

## Read first

- `docs/decisions/0027-outlet-void-and-discharge-rate.md` (the outlet fix and every
  measurement after it), `0026` (dissipation study before the fix), `0023` (shear cap)
- `docs/wp/WP-3.3q.md` log; `docs/results/silo-prototype.md`
- `prototype/void2d.py` (rules 8, 10, 10b, 11 and the `DISS_N`, `DISS_T` knobs)

## Build

- First find what sets the E dependence. It is not a coarse-mesh dissipation of order
  rho c h: at E = 1 MPa the rate falls with refinement (73 kg/s/m at 10 mm, 61 at 5 mm)
  while at 0.1 MPa it does not move (81, 82). Candidates: elastic dilation of the soft
  bed near the outlet (a physical effect that stiffer material lacks), the number of
  steps (dt scales with h / c), the void-first closure (rule 4) in the stretching
  outlet zone. Measure density and al in the outlet zone against E, and E = 10 MPa at
  5 mm (10 mm: 54 kg/s/m).
- Then a damping of the elastic and odd-even modes that does not resist plastic flow
  (for example jump terms reduced only at faces where both cells are at yield, or an
  artificial viscosity on the volume change only). Measure each with
  `discharge2d.py --run --cell 0.01 --E 1e7 --t-end 0.6` (about 15 minutes): rate,
  0.1 s window variation, largest speed, corner speed (averaged, `averaged_pattern`).
- If the stiff rate still depends on the numerics, mu(I) friction (Jop et al. 2006) in
  `constitutive/batch.py` as an option, and report both.
- Keep the void, collapse, settle, discharge and rate tests passing unchanged.

- Also explain the fine-mesh rise at E = 0.1 MPa (D = 0.06 m: 81.4 / 82.2 / 93.3
  kg/s/m at 10 / 5 / 2.5 mm).

## Done when

A new `tests/test_proto_silo_stiffness.sh` (10 mm, D = 0.06 m) passes: rate at
E = 0.1, 1 and 10 MPa within 10 percent of each other, 0.1 s windows within 10 percent,
averaged corner speed below 1 percent of the outlet speed at each E. Reported: 5 against
2.5 mm at E = 0.1 MPa within 10 percent.

## Out of scope

- Cohesion (3.5), C++ (3.3, 3.4), the outlet (done, rule 10b).

## Log

(none yet)
