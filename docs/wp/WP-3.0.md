# WP 3.0: Silo benchmark definition and reference values

State: not started
Branch: `wp-3.0-silo-reference`

## Goal

`docs/spec/silo-reference.md`: the numbers and curves the silo simulations (WP 3.1p to
3.5) are judged against, each with a source that was actually opened and read.

## Read first

- `docs/spec/target-problems.md` (T0), `docs/decisions/0016-silo-discharge-first.md`

## Build

- Janssen's theory for a 2D (plane) bin: vertical and wall stress against depth for the
  T0 geometry and material, with the stress ratio K (state which K and why: Jaky,
  active Rankine, or measured ranges) and the wall friction; saturation stress and the
  depth at which 90 percent of it is reached. A small script
  `tools/silo_reference.py` printing these for any inputs.
- Discharge rate: the Beverloo/Hagen law in its 2D (slot) form, the range of its
  constant reported for continuum simulations (which have no particle diameter), and
  the expected rate in kg/s per metre depth for D = 0.04, 0.06, 0.08 m. The condition
  for rate independent of head.
- Published continuum simulations of flat-bottomed silo discharge (finite element,
  MPM, mu(I) rheology, hypoplastic; for example Staron, Lagree and Popinet's "granular
  silo as a continuum plastic flow", and Mohr-Coulomb or Drucker-Prager Eulerian FE
  studies): what rate constant, flow pattern and stagnant-zone angle they found, and
  with which constitutive law. Note where Mohr-Coulomb with perfect plasticity is known
  to differ from experiments.
- Granular column collapse in 2D (for WP 3.2p): run-out and final height against the
  initial aspect ratio from experiments (Lube et al., Lajeunesse et al.) with the
  fitted laws and their ranges.
- For each number: source, link, and how tight a simulation should be expected to match
  (a tolerance to use in the tests of 3.1p, 3.2p, 3.3p and 3.5).

## Done when

- `docs/spec/silo-reference.md` exists with every value sourced (title, authors, year,
  link, and "opened on" date); anything not found in an opened source is marked
  "not verified".
- `tools/silo_reference.py` runs in CI inside `tests/test_silo_reference.sh` and
  reproduces the tabulated Janssen and Beverloo numbers.
- `docs/STATUS.md` and the briefs of 3.1p, 3.2p, 3.3p updated where a tolerance changes.

## Out of scope

- Any solver code.

## Log

(none yet)
