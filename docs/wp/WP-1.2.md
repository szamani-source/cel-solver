# WP 1.2: 1D plasticity and history remap

State: done, waiting for review
Branch: `claude/project-thread-a7fstf` (the session's assigned branch; the planned name was `wp-1.2-1d-plasticity`)

## Goal

The 1D prototype carries an elastic-plastic material: the stress update calls the WP 1.4
models, and plastic strain (a history variable) is remapped with the material. This is
the first test of the hardest CEL idea: history must move with the material, not smear.

## Read first

- `docs/decisions/0005-cell-centred-lagrange-step.md` (the Lagrange step and remap in use)
- `docs/decisions/0004-constitutive-prototype-conventions.md` (model signature, Voigt order)
- `prototype/lr1d.py` (class `Bar1D`), `prototype/constitutive/models.py` (`von_mises`,
  `drucker_prager`)
- `docs/wp/WP-1.1.md` log, for what was tried and why

## Build

- Extend `prototype/lr1d.py`, or add `prototype/lr1d_plastic.py` that reuses it, so that:
  - the bar is in uniaxial strain (lateral strain zero), with the full 6-component stress
    kept per cell and updated by a WP 1.4 model from the axial strain increment;
  - the wave speeds used by the Riemann solver and the time step are the elastic ones;
  - history (`eps_p`, `eps_p_eq`) and all stress components are remapped conservatively
    (mass-weighted) with the same limiter as the other fields.
- `tests/test_proto_1d_plastic.sh`
- A decision note if the Godunov step needs a change for the two-wave (elastic precursor,
  plastic wave) structure.

## Done when

`tests/test_proto_1d_plastic.sh` passes, checking with von Mises (no hardening):

- Step loading above the Hugoniot elastic limit gives two waves. Elastic precursor speed
  within 1 percent of sqrt((K + 4G/3)/rho); plastic wave speed within 2 percent of
  sqrt(K/rho); precursor amplitude within 2 percent of the Hugoniot elastic limit
  sigma_HEL = Y (K + 4G/3) / (2G).
- No new extrema above 1 percent of the jump in stress or velocity.
- Streaming test: a patch of pre-strained material (non-zero `eps_p_eq`) carried at
  0.25 c across at least 50 cells with zero loading keeps its total (mass-weighted)
  plastic strain to round-off, creates no plastic strain elsewhere, and its edges stay
  within 4 cells wide.
- Mass and momentum conserved to round-off. Existing tests still pass unchanged.

Also run, and record in the log (not a pass criterion): the same step loading with
Drucker-Prager, to see that a pressure-dependent model runs stably.

## Out of scope

- 2D (WP 1.3), more than one material or void (Phase 3), hardening laws beyond what
  WP 1.4 has, the extremum-preserving limiter (known issue, only if this package needs it).

## Log

### 2026-11-14

Done:
- `prototype/lr1d_plastic.py`: class `PlasticBar1D(Bar1D)`. Uniaxial strain, 6-component
  stress and history per cell, any WP 1.4 model by name. `lr1d.py` got one hook,
  `Bar1D.update_stress`; its own output is bit-identical to before.
- `tests/test_proto_1d_plastic.sh` (runs `lr1d_plastic.py --verify`, about 25 s). All
  checks pass at 1000 cells, 0.25 c_e streaming: precursor speed error 0.04 percent,
  plastic wave 0.09 percent, precursor amplitude 0.18 percent (15 plateau cells); extrema
  after start-up 0.13 percent stress, 0.11 percent velocity; mass and momentum change
  below 2e-16; streaming patch keeps total plastic strain to 3e-14 after 120 cells, edges
  2 cells (10-90 percent), 5 cells (1-99 percent), no growth from 60 to 120 cells.
- Decision 0008: elastic impedance kept, superbee for history, start-up overshoot.
- Drucker-Prager run (`lr1d_plastic.py --dp`, not a pass criterion): E 20 MPa, nu 0.2,
  phi 20, c 5 kPa, psi 0. Stable. Precursor speed 105.43 (c_e 105.41), plastic wave
  91.86 m/s (model probe 91.73), precursor 28.9 kPa vs probe HEL 29.1 kPa (0.9 percent
  low); extrema after start-up 0.34 percent, start-up 2.5 percent; tension side sits at
  the apex (13.74 kPa = xi c/eta) with no blow-up; mass and momentum to round-off.

Learned:
- With nu 0.3 and phi 30, Drucker-Prager never yields in uniaxial-strain compression
  (elastic K0 0.43 is inside the cone). The --dp run uses nu 0.2, phi 20.
- An impulsive start makes a 2.7 percent overshoot at the impact plane for the first
  5 to 10 cells of wave travel. The elastic bar does it too (1.7 percent), so it is the
  WP 1.1 scheme under impulsive loading. Inside moving waves extrema stay near 0.1 percent.
- MC remap of a discontinuous history field keeps widening slowly; superbee does not.
- `uniaxial_response` probes any model for its HEL and plastic wave speed; for von Mises
  it reproduces the closed forms to 1e-9.

Update: Sia chose to check extrema after the start-up (decision 0008 item 3, accepted).

Next: Sia reviews and merges PR #6. Then WP 1.3.
