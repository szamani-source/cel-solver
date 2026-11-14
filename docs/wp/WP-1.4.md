# WP 1.4: Single-point constitutive driver (elastic, J2, Drucker-Prager, Mohr-Coulomb)

State: not started
Branch: `wp-1.4-constitutive-driver`

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
- `tests/test_proto_constitutive.sh`

## Done when

- Triaxial compression and extension strengths match Mohr-Coulomb closed form within 0.1 percent.
- Drucker-Prager fitted to the compression meridian gives the same triaxial compression strength.
- Returned stress lies on or inside the yield surface for 10,000 random strain increments,
  including apex and corner cases.
- Volume change in simple shear matches the dilation angle.

## Out of scope

- C++ port (WP 2.2, 2.5a, 2.5b), softening, rate effects.

## Log

(none yet)
