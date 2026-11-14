# WP 1.4: Single-point constitutive driver (elastic, J2, Drucker-Prager, Mohr-Coulomb)

State: done (CI green on PR #2, waiting for merge)
Branch: `claude/project-thread-w8f1nv` (the session's assigned branch; the brief had proposed `wp-1.4-constitutive-driver`)

## Goal

A Python single-point driver (`prototype/constitutive/`) with the four models in
`docs/constitutive-models.md`, each verified against closed-form results. These become the
reference the C++ library is tested against.

## Read first

- `docs/constitutive-models.md`
- de Souza Neto, Peric and Owen, "Computational Methods for Plasticity", chapters on
  Mohr-Coulomb and Drucker-Prager return mapping

## Build

- `prototype/constitutive/models.py`: one function per model, same signature
  (strain increment, stress, history, parameters) -> (stress, history).
- `prototype/constitutive/driver.py`: strain-driven and mixed-control paths
  (uniaxial strain, triaxial compression and extension, simple shear).
- `prototype/constitutive/check_models.py`: the checks, run by
  `tests/test_proto_constitutive.sh`

## Done when

- Triaxial compression and extension strengths match Mohr-Coulomb closed form within 0.1 percent.
- Drucker-Prager fitted to the compression meridian gives the same triaxial compression strength.
- Returned stress lies on or inside the yield surface for 10,000 random strain increments,
  including apex and corner cases.
- Volume change in simple shear matches the dilation angle.

## Out of scope

- C++ port (WP 2.2, 2.5a, 2.5b), softening, rate effects.

## Log

### 2026-11-14

Done:
- Wrote `tests/test_proto_constitutive.sh` and `check_models.py` first, then `models.py`
  (linear elastic, J2 with linear hardening, Drucker-Prager with compression, extension
  and plane-strain fits, Mohr-Coulomb with main plane, corner and apex returns) and
  `driver.py` (mixed control by Newton with a finite-difference Jacobian).
- All 20 checks pass locally (about 7 s, Python 3 and numpy only):
  MC triaxial compression and extension strengths match closed form to about 1e-11;
  DP compression fit gives the MC compression strength and the extension fit the
  extension strength; 10,000 random increments per case (6 cases: MC, MC c = 0,
  MC associated, DP compression, DP extension c = 0, DP plane strain) end on or inside
  the surface (worst f/scale 4e-14), with more than 1,500 apex returns per cohesionless
  case; simple-shear dilatancy d(eps_v)/d(gamma) = tan(psi) within 0.01 percent for MC
  and the plane-strain-fitted DP.
- Mutation check: swapping the two MC corner branches makes the triaxial checks fail,
  so the corner logic is tested.
- Conventions and algorithm choices recorded in decision 0004.

Learned:
- In plane-strain simple shear at steady plastic flow MC gives d(eps_v)/d(gamma) =
  tan(psi) exactly; DP gives eta_bar / sqrt(1 - 4 eta_bar^2 / 3), which equals tan(psi)
  only for the plane-strain fit.
- Stale `__pycache__` can hide an edit of the same size made within one second; it is
  now ignored by git.

CI: green on PR #2 (build-and-test, 2026-11-14). No separate tension cut-off yet
(see decision 0004).

Next: WP 2.2 / 2.5a / 2.5b port these models to C++ and reuse `check_models.py` values
as reference numbers.
