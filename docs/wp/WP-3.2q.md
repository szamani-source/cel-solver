# WP 3.2q: Prototype column collapse (split from WP 3.2p)

State: not started
Branch: thread branch

## Goal

The 2D granular column collapse in the void prototype passes the decision 0019 checks, or
the reason it cannot is shown with numbers and the bands are revisited in a decision
note (never silently).

## Read first

- `docs/decisions/0021-void-treatment.md` (section "Column collapse"), `0019`, `0020`
- `docs/wp/WP-3.2p.md` log; `prototype/void2d.py` (`column_collapse`, `collapse_check`,
  the `DISS` study knob on the Riemann jump terms)
- `docs/spec/silo-reference.md` section 4

## Build

- Find what sets the run-out: separate the numerical dissipation of the acoustic Riemann
  solver (run-out depends on E, which a rate-independent model should not) from the
  physics of phi = 35 deg. Suggested order: an E series (1, 0.3, 0.1, 0.03 MPa) at 10
  and 20 cells; a low-dissipation flux for slow flow (for example the jump terms scaled
  with a local speed ratio, or an artificial-viscosity form in place of the full
  acoustic upwinding), checked against the WP 1.3a/2.3 wave tests so they stay passing.
- The blunt front (2 to 4 cells high, 50 to 60 deg per cell): try second-order
  reconstruction next to partly filled cells (a plain "cell and neighbours hold
  material" mask blew up the falling block; partly filled cells need their own slopes).
- If phi = 35 deg truly runs out shorter than 0.5 times 1.2 a, show it (for example
  an energy or wedge estimate, or a published rate-independent Mohr-Coulomb result) and
  write the revised band as a decision.
- Then the 3.3p discharge rate should be checked for the same E dependence.

## Done when

`tests/manual/test_proto_collapse.sh` passes and is moved to `tests/` (so CI runs it):
a = 0.5 keeps its top within 1 cell; run-out 0.30 to 0.90 (a = 0.5) and 1.2 to 3.6
(a = 2); no deposit slope steeper than 37 deg; mass to 1e-10. Existing tests unchanged.

## Out of scope

- C++; the outlet (3.3p); cohesion.

## Log

(none yet)
