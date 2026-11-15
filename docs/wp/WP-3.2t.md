# WP 3.2t: One rule set for the collapse and the discharge (split from WP 3.3r)

State: closed without meeting its done-when (2026-11-15, decision 0029)
Branch: thread branch

## Goal

Rules 12 (normal jump cap) and 13 (dispersed interior void) of `prototype/void2d.py`
are either on everywhere or replaced by something that is, so the column collapse and
the silo discharge use the same numerics, which the C++ Lagrange step can then copy.

## Read first

- `docs/decisions/0028-normal-jump-cap-and-dispersed-void.md`, `0023` (shear cap),
  `0025` (collapse toe, creep)
- `docs/wp/WP-3.3r.md` log
- `prototype/void2d.py` rules 4, 8, 12, 13; `column_collapse`

## Build

- With both rules on, the a = 0.5 collapse (10 cells) runs out 0.8 / 1.05 / 1.2 at
  E = 0.1 / 1 / 10 MPa (band 0.30 to 0.90; off: 0.45 / 0.40 / 0.25). Find which rule
  does it (run each alone) and why the collapse does not stop: the jump terms were the
  only damping of the slow front (decision 0025 found the a = 0.5 front creeping at
  faces held at yield).
- Candidates: rule 12 only where the cell is flowing fast (shear rate above a small
  threshold), a stick condition for slow faces (`VSTICK`, decision 0025), or mu(I)
  (rate-dependent friction, Jop et al. 2006) which also stops slow creep.
- Keep every existing test passing unchanged.

## Done when

`tests/test_proto_collapse.sh` passes with rules 12 and 13 on (or their replacement),
with run-outs at E = 1 and 10 MPa within 0.1 of each other, and
`tests/test_proto_silo_stiffness.sh` still passes.

## Out of scope

- C++ (WP 3.1 to 3.4), cohesion (3.5).

## Log

### 2026-11-15

- Each rule alone breaks the collapse: rule 13 only 0.80 / 0.85 / runs away at
  0.1 / 1 / 10 MPa; rule 12 only 0.95 / runs away (a = 0.5, 10 cells).
- Rule 13: void spread by the remap diffuses through the bulk as stress-free porosity
  (90 interior partly filled cells at 1 MPa against under 7). Keeping isochoric stretch
  in the material (`STRETCH_ISO`) cures it (0.40 at 1 MPa) but slows the discharge
  (76 against 87 kg/s/m at 1 MPa).
- Rule 12: the plain cap turns the column hydrostatic (deviator ratchets down at capped
  faces of elastic cells). Capping only at yield fixes that and gives 0.35 at 1 MPa,
  but the 10 MPa column never stops (time-step dependent: the jump is what keeps the
  split step stable at faces at yield), and the discharge falls to 48 to 52 kg/s/m at
  1 MPa. Viscous floor, rate threshold, full cells only, volumetric part kept,
  geometric strain and an odd-even jump all fail one case or the other (tables in
  decision 0029).
- So the discharge's E independence (0028) leans on capping faces whose cells are
  still elastic. Defaults unchanged; study knobs kept, all off.
- `tests/test_proto_collapse.sh`, `test_proto_void.sh`, `test_proto_silo_discharge.sh`
  pass unchanged with the knobs in place.
- Next: the C++ silo packages copy the full jump terms (no rules 12 or 13), then measure
  the discharge rate against E and the collapse with the same code (decision 0029).
