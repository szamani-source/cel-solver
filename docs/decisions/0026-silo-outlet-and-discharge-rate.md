# 0026 Silo outlet boundary; discharge rate set by Riemann dissipation (split to WP 3.3q)

Date: 2026-11-15. Status: accepted for the prototype (WP 3.3p).

## Context

WP 3.3p runs target T0 end to end in `prototype/discharge2d.py` on `void2d.VoidBox`:
tank 0.40 m wide, fill 0.60 m, walls and floor with 20 deg friction, Mohr-Coulomb c = 0,
phi 35 deg, psi 0. Settling as WP 3.1p (gravity ramp 0.1 s, damping 250 1/s, decision
0020), for a fixed 0.3 s; then the outlet opens and damping is switched off. The brief
asked to measure the rate's dependence on E first (decisions 0021, 0023).

## Decision 1: the outlet (rule 10 of `void2d.py`)

`void_sweep(..., outlet=mask)` marks rows (columns of the tank) whose bottom face is an
outflow face. Those rows get, in the gravity sweep only:

- empty ghosts for the volume fraction (material leans towards the tank, not out);
- free-end ghost parities and the free boundary Riemann solution: zero normal and
  shear traction (s_hat = -G at the face, which is 0 at the bottom);
- a remap volume flux through the face equal to the material swept across it (the
  same slab overlap as interior faces), at the donor's Lagrangian density. That mass
  leaves the domain and is added to `VoidBox.discharged` (averaged over the two sweep
  orderings, as the state is).

The outlet opens instantly (`VoidBox.outlet = mask`); no start-up spike was seen, so the
ramp option of the brief was not built. Below the outlet nothing is modelled: the
material leaves at the face. Result (10 mm, E = 0.1 MPa): nothing passes the closed
outlet (exactly 0), and mass in tank plus discharged mass is constant to 1e-14.

## Decision 2: the rate is not a result yet; the rate checks move to WP 3.3q

The steady rate is set by the numerical dissipation of the acoustic Riemann solver, not
by the material. Measurements, D = 0.06 m, 10 mm cells, rate over 0.05 to 0.4 s after
opening (Beverloo with C = 1.45 gives 120 kg/s/m):

| E (MPa) | 0.01 | 0.1 | 1 | 10 (T0) |
| --- | --- | --- | --- | --- |
| Rate (kg/s/m) | 47 rising to 274 (not steady) | 38 | 17 | 6 |

The jump terms z (vR - vL) / 2 act as a viscosity of order rho c h |grad v| in the
converging flow above the outlet; it grows with sqrt(E) and, unlike the one-cell shear
bands of decision 0023, it does fall with the cell size (E = 0.1 MPa, rate over 0.2 to
0.5 s: 36.8 kg/s/m at 5 mm and 45.2 at 2.5 mm, 23 percent apart; see
`docs/results/silo-prototype.md`). The free-surface rule matters too: WP 3.2r's rule 9
(merged during this package) lowered every rate by 20 to 40 percent. The E table above
was measured before that merge.

Tried and rejected (E = 1 MPa, 10 mm, unless stated):

| Change | Rate (kg/s/m) | Kinetic / strain energy | Largest speed (m/s) |
| --- | --- | --- | --- |
| None | 17 | 0.04 to 0.17 | 0.45 |
| Shear cap (0023) off | 9 | | |
| Both jump terms times 0.3 | 27 | | |
| Both jump terms times 0.1 | 58 to 78 | 1.4 to 2.5 | 4 to 7 |
| Normal jump term times the face Mach number min(1, abs(v)/c) | 30 | | |
| Both jump terms times the Mach number (Thornber et al. 2008 low-Mach fix) | 42 to 55 | 2.6 to 6.6 | 3.5 to 9.5 |
| Same with a floor of 0.03 / 0.1 on the factor | 36 to 63 / 24 to 49 | 2 to 6 / 1 to 3.7 | 4 to 10 |

Every reduction of the jump terms raises the rate but lets the flow turn noisy: speeds
of 4 to 10 m/s (free fall from the fill height is 3.4 m/s), kinetic energy above the
elastic strain energy, corner speeds 3 to 90 percent of the outlet speed. With the
Mach scaling during settling the column rang (kinetic energy 10 times the strain
energy at 0.3 s). The Mach-scaled rate also did not converge in E (40, 48, 58, 105, 165
kg/s/m for E = 0.1, 1, 3, 10, 30 MPa). The jump terms are the only damping of the
elastic and odd-even modes (decision 0023), so they cannot simply be removed. The Mach
scaling is kept in the code as `low_mach` (rule 11, off by default) for WP 3.3q.

A second, physical reason is recorded, not tuned: rate-independent Coulomb friction
does not select a unique steady flow (`docs/spec/silo-reference.md` section 3, S4; the
continuum values of C all come from mu(I) rheology, where the rate dependence selects
it). So even a dissipation-free scheme may leave the rate set by the numerics.

## Consequences

- `tests/test_proto_silo_discharge.sh` (CI) checks what is done: the outlet mechanics
  (no leak, mass balance 1e-8, outflow, head independence within 10 percent: 0.975).
- The WP 3.3p rate checks (steady rate, exponent, C, corners, crater) are in
  `tests/manual/test_proto_silo_rate.sh` and fail: at E = 0.1 MPa, 10 mm, exponent 2.2
  (band 1.3 to 1.7), C 0.38 / 0.51 / 0.61 for D = 0.04 / 0.06 / 0.08 m (band 1.0 to 1.9;
  before the WP 3.2r merge); at 5 mm with rule 9, exponent 2.26 and C 0.32 / 0.44 /
  0.54. The rate is steady (within 7 percent at 5 mm), the corners are stagnant and a
  crater forms: the flow pattern is right, the magnitude is a third to a half of
  Beverloo's.
  They are the done-when of WP 3.3q (`docs/wp/WP-3.3q.md`).
- WP 3.3q options, in order: (1) a damping that acts on elastic and odd-even modes but
  not on the plastic flow (for example jump terms scaled only where both cells are at
  yield, or a separate artificial viscosity on the velocity divergence); (2) mu(I)
  friction (Jop et al. 2006) in `constitutive/batch.py`, which gives the flow a rate
  dependence that selects the steady rate and is what every continuum reference used;
  (3) both. Use E = 0.1 MPa while developing: it is 10 times cheaper than T0's E.
- The C++ outlet (WP 3.3) can copy rule 10 as it is.
