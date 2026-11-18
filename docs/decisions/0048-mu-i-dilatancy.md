# 0048 Dilatancy in the mu(I) law (WP 3.6g)

Date: 2026-11-18. Status: accepted (option kept, off by default; done-when not met).

## Context

Decision 0046 traced the 1 MPa dip in the T0 discharge rate (bands 65.8 / 52.0 / 60.8
kg/s/m at 0.1 / 1 / 10 MPa, 10 mm) to how loosened cells above the outlet are treated
(full or partly filled), and proposed making the loosening a state of the material:
a compressible mu(I)/phi(I) law. Sia approved the change ("proceed", 2026-11-18).
Done-when: the three bands within 10 % of each other at 10 mm, with
`tests/manual/test_inclined_layer.sh` still passing.

## What was built

Optional entries `phi0`, `phiC`, `bPhi`, `KDil` (all four or none; mu(I) needed) in
`mohrCoulomb` and `druckerPrager`. The dilation angle follows the state instead of the
fixed `psi` (Roux and Radjai 1998; Pailha and Pouliquen 2009):

    tan(psi) = KDil (phiS - phiEq(I)),  phiEq(I) = phiC - bPhi I,
    phiS = phi0 exp(-tr(epsP))

- phiS is the solid fraction from the plastic volumetric strain already in the history
  (no new field); phi0 is its value at the fill (rho / rhoS).
- In the return at friction mu, I is the Jop law's I(mu), so psi and mu agree at the
  root of the existing Illinois solve; with dilatancy the mu_s return also uses the
  state psi. |psi| <= phi_m, and contraction is bounded so the plastic multiplier
  stays positive (MC: 4G + 4 sin(phi) sin(psi) (lambda + G) >= 2G; DP:
  G + K eta etaBar >= G/2).
- The nearly-full snap (decisions 0044, 0045) compares with the dilated density
  rho0 exp(-tr(epsP)) (1 + p/K).
- Without the entries the models are unchanged bit for bit (`tests/test_mu_i.sh`,
  `tests/test_constitutive_cpp.sh` pass).

Kernel check `tests/test_dilatancy.sh` (CI, 2 s): plane-strain shear box at 5 kPa,
dense (phi0 0.679) and loose (0.55) starts both reach phiEq(I) to 1e-9; steady friction
equals plain mu(I) (checked to 1e-6, equal to 8 digits); first plastic step d(eps_v^p)/d(gamma^p) = 0.308, inside
the sin(psi) bracket; Drucker-Prager reaches phiEq too.

## Results

Placeholders: bPhi 0.30 (da Cruz et al. 2005: 0.31), KDil 4 (Pailha and Pouliquen
2009: 4.09), phi0 = 1800/2650 = 0.679. Tank bands, 10 mm, three runs each
(`tests/manual/test_discharge_scatter.sh`, entries through `CEL_EXTRA_MATERIAL`):

| Variant | 0.1 MPa | 1 MPa | 10 MPa |
| --- | --- | --- | --- |
| Without dilatancy (decision 0046) | 65.8 (65.5 to 66.3) | 52.0 (50.6 to 53.8) | 60.8 (56.5 to 67.7) |
| Dense fill, phiC 0.60 (psi about 18 deg at the start) | 53.6 (52.8 to 54.1) | 35.9 (33.1 to 40.8) | 32.1 (29.3 to 34.3) |
| Fill at the critical state, phiC 0.679 (dilation only from bPhi I) | 61.5 (61.4 to 61.6) | 42.3 (41.6 to 43.2) | 45.8 (43.9 to 47.1) |
| Same, KDil 1e-9 (control: only the snap change) | | 50.4 (49.0 to 51.8) | |

Inclined layer (5 mm, started from the exact steady state, phiC set so that
phiEq(I_theta = 0.0683) = phi0, i.e. psi = 0 in the exact solution): profile errors
4.6 / 2.6 / 6.7 % at 0.1 / 1 / 10 MPa (limit 5 %; without dilatancy 3.6 / 0.9 / 0.8 %).

## Findings

1. Dilatancy slows the tank at every E and widens the gap: the 1 MPa band moves from
   21 % below 0.1 MPa to 31 % below, and 10 MPa falls with it. The control run puts
   all of it on the dilation, not on the snap change.
2. Mechanism (reasoned, not measured cell by cell): dilation in the flowing zone is
   confined by the material around it, so it raises the pressure, and with it the
   friction, by about K times the dilation strain until the extra volume has moved
   away. At 0.1 MPa the elastic strains (about 5 %) swallow a dilation of a few
   percent; at 1 and 10 MPa they do not. So in this explicit scheme dilatancy adds an
   E-dependent resistance; it does not remove one.
3. The loosening that sets the rate (decision 0046) is opened at zero stress, in void,
   by material pulled apart above the outlet. Plastic dilatancy acts only while the
   material carries stress, so it does not reach that loosening.
4. The layer is no longer uniform in I near the floor and the top cell, so those
   cells dilate or contract and the profile drifts by a few percent.

## Decision

- Keep the dilatancy entries as an option, off by default and in no tutorial. It is a
  correct material-point model (checked) and wanted later for the collapse shear bands
  (decision 0038), where the material is confined less.
- Close WP 3.6g with its done-when not met. The 1 MPa dip stays a known issue
  (decision 0046 rule unchanged).
- Proposed next (not started): treat the full / partly filled switch as a smooth
  blend over void fractions 1e-4 to 1e-2 instead of a threshold, so the silo rate
  stops depending on which side of the line loosened cells fall (decision 0046,
  finding 3). That is a numerical change, not a model change.

## Consequences

- New optional material entries `phi0`, `phiC`, `bPhi`, `KDil` (documented in
  `celMuI.H` and here; no tutorial uses them yet).
- Placeholders are not Sia's material. phiC for sand is about 0.57 to 0.62; T0's fill
  density (1800 kg/m3 with rhoS 2650, phi0 0.68) is denser than most sands' densest
  packing, worth checking with Sia if dilatancy is used.
