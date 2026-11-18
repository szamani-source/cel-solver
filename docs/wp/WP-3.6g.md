# WP 3.6g: dilatancy in the mu(I) law

State: closed, done-when not met (2026-11-18)
Branch: `claude/project-thread-4zdz5u`

## Goal

Make loosening a state of the material (compressible mu(I)/phi(I)) so the T0 tank rate
agrees across E, without breaking the inclined layer (decision 0046 proposal).

## Read first

Decisions 0037, 0044, 0045, 0046, 0048; `src/celConstitutive/kernel/celMuI.H`.

## Done when

`tests/manual/test_discharge_scatter.sh` with dilatancy: 1 and 10 MPa bands within
10 % of 0.1 MPa at 10 mm, and `tests/manual/test_inclined_layer.sh` still passes.
Kernel: `tests/test_dilatancy.sh` (CI).

## Log

- 2026-11-18: One session. Built the state-dependent dilation angle (entries phi0,
  phiC, bPhi, KDil) in MC and DP; kernel check passes. Tank bands: dense fill 53.6 /
  35.9 / 32.1, fill at critical 61.5 / 42.3 / 45.8 kg/s/m (without 65.8 / 52.0 / 60.8):
  slower and more E dependent. Layer 4.6 / 2.6 / 6.7 % (fails at 0.1 and 10 MPa).
  Confined dilation raises the pressure by about K times the dilation strain, so it
  adds an E-dependent resistance; the loosening that sets the rate happens at zero
  stress, which plastic dilatancy does not reach. Option kept, off by default.
  Next step proposed: smooth blend of the full / partly filled rules (decision 0048).
