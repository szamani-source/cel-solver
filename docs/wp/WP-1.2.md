# WP 1.2: 1D plasticity and history remap

State: not started
Branch: `wp-1.2-1d-plasticity`

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

(none yet)
