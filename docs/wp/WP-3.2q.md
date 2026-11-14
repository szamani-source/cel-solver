# WP 3.2q: Prototype column collapse (split from WP 3.2p)

State: partly done (run-out passes; slope split to WP 3.2r), 2026-11-14
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

### 2026-11-14

- E series (a = 0.5, 10 cells): 0.50, 0.40, 0.35, 0.25, 0.25, 0.20 for E = 0.01, 0.03,
  0.1, 1, 3, 10 MPa. Not the time step (smaller CFL gives longer run-out), not the
  reconstruction order next to the surface, mesh independent (10 and 20 cells agree).
  Scaling the Riemann jump terms by 0.3 gives 0.40: the jump term z dv / 2 in a
  one-cell shear band is a mesh-independent viscosity growing with sqrt(E).
- Fix (decision 0023, rule 8, `shear_cap`): the shear jump term is cut back so the face
  stress stays inside Mohr-Coulomb; capped faces slip like a frictional wall. a = 0.5
  0.40 (10 and 20 cells), top 0.06 cell lower; a = 2 1.6. Both run-out bands pass;
  `tests/test_proto_void.sh` unchanged and passing.
- Rejected: Coulomb limit on the face plane only (no effect); capping the normal jump
  term (any form) or the whole shear traction (unstable or the top sinks).
- Not passed: deposit slope 45 deg (a = 0.5) and 59 deg (a = 2) per cell, two-cell 35
  and 53. Cause: the 1D slab rule (0021 rule 2) keeps thin top layers from moving
  sideways, so column heights snap to whole cells. A 2D-normal slab rule smooths it
  without the cap (0.35) but blew up with the cap. Split to WP 3.2r; the test stays in
  `tests/manual/`. Residual E dependence (0.45, 0.40, 0.25 at 0.1, 1, 10 MPa) from the
  uncapped normal jump term is recorded in 0023 for 3.3p.
