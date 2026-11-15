# WP 3.3q: Prototype silo discharge rate (split from WP 3.3p)

State: not started (split 2026-11-15, decision 0026)
Branch: thread branch

## Goal

The T0 discharge in `prototype/discharge2d.py` gives a steady rate that follows the 2D
Beverloo law (exponent 1.5 +- 0.2, C in 1.0 to 1.9) and does not depend on E or on the
mesh, with no spurious noise in the stagnant zones.

## Read first

- `docs/decisions/0026-silo-outlet-and-discharge-rate.md` (all measurements so far),
  `0023` (shear cap), `0021` (void rules)
- `docs/spec/silo-reference.md` section 3 (C, mu(I) caveat), section 6 (tolerances)
- `docs/results/silo-prototype.md`; `prototype/discharge2d.py`, `prototype/void2d.py`
  (rules 8 to 11)

## Build

- Find a damping for the elastic and odd-even modes that does not act as a viscosity on
  the plastic flow. Candidates: jump terms scaled by the Mach number only at faces where
  both cells are at yield (rule 11 is the all-faces version, which turned noisy); an
  artificial viscosity on the velocity divergence; something from the literature on
  low-dissipation Godunov schemes for solids. Measure each with
  `discharge2d.py --run --cell 0.01 --E 1e6 --t-end 0.3` and record rate, kinetic over
  strain energy, largest speed and corner speed (as in the decision 0026 table).
- If the rate still depends on the numerics once the noise is gone, add mu(I) friction
  (Jop et al. 2006; parameters from S2 or S4 in `silo-reference.md`) as an option of the
  Mohr-Coulomb model in `constitutive/batch.py`, and report both.
- Keep `tests/test_proto_void.sh`, `tests/test_proto_silo_settle.sh` and
  `tests/test_proto_silo_discharge.sh` passing unchanged.

## Done when

`tests/manual/test_proto_silo_rate.sh` passes and moves to `tests/` (10 mm cells,
E = 0.1 MPa): steady rate (0.1 s windows within 10 percent), exponent 1.3 to 1.7, C 1.0
to 1.9, corner speed below 1 percent of the outlet speed, a crater, head independence 10
percent. Plus, reported in `docs/results/silo-prototype.md`: rate at E = 0.1, 1, 10 MPa
within 10 percent, and 5 mm against 2.5 mm within 10 percent.

## Out of scope

- Cohesion (3.5), C++ (3.3, 3.4), the collapse deposit slope (3.2r).

## Log

(none yet)
