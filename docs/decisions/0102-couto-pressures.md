# 0102: Couto silo pressures: T0 rule set for both models, results

Date: 2026-11-19 (WP H.3)

## Context

Case A of the conical hopper validation (decision 0100) compares the wall pressure of
the Couto et al. (2013) wheat silo, at rest and in discharge, with the experiment.
P1 (Zheng and Yu) used plain, rate-independent Mohr-Coulomb in Abaqus.

## Findings

- **Plain MC and DP (no mu(I)) are not usable once the outlet opens.** Plain MC flows,
  and its pressures averaged over time land near the experiment. But single cells swing
  from 0 to 7 kPa, and stresses reach 0.9 MPa with speeds of 18 m/s, against a physical
  25 kPa and about 2 m/s. Plain DP blows up: 5 MPa and 40 m/s within 1 s of discharge.
  This is the known need for the T0 rule set on free-surface tanks (WP 3.8b, 4.3b).
- **T0 rule set (mu(I) + flowImpedance + normalCap + shearCap + spreadLayers).**
  `flowImpedance` was allowed for Mohr-Coulomb only. It now accepts druckerPrager too
  (`eulerianSolid.C`): the face rules cap jumps at the Mohr-Coulomb yield of phi for
  both models. The mu(I) term uses Jop's dmu 0.26 and I0 0.279 (as T0), with wheat
  d = 4 mm and rhoS = 1380 kg/m^3 (856 / 0.62). At the silo's I (below about 0.03),
  it raises the friction by less than 0.03, about 4 % of tan(phi). P1 had no rate term.
- With the rule set MC runs clean (stress below 40 kPa, speed about 1 m/s). DP mostly
  runs too, but it is fragile. Runs that differ only in round-off went up to 680 kPa
  and 8 m/s in one case and stayed at 11 m/s in another. In every case the trouble
  starts at partly filled cells near the outlet or at the hopper's cut cells. No
  tension was found there (largest principal stress 0), so it is not the missing
  DP tension cut-off. The DP fit (compression, extension, plane strain) does not cure it.

## Results (manual `tests/manual/test_hopper_pressure.sh`, 25 mm, 12 min per model)

| Cell | At rest MC / DP | Experiment at rest | Discharge mean MC / DP | Experiment discharge |
| --- | --- | --- | --- | --- |
| GA7 | 1.36 / 1.36 | 1.39-1.69 | 1.92 / 1.25 | 1.43-1.74 |
| GA6 | 2.12 / 2.12 | 2.47-2.66 | 2.60 / 1.95 | 3.35-3.63 |
| GA5 | 2.82 / 2.82 | 2.78-3.24 | 3.27 / 2.76 | 2.88-3.72 |
| GA4 | 3.48 / 3.48 | 3.17-3.62 | 3.78 / 3.61 | 3.21-4.01 |
| GA3 | 4.04 / 4.04 | 3.67-4.16 | 4.22 / 3.85 | 4.42-5.59 |
| GA2 | 4.53 / 4.52 | 4.07-4.46 | 5.13 / 4.10 | 5.81-6.45 |
| GA1 | 4.79 / 4.80 | 4.00-4.49 | 6.95 / 6.61 | 7.14-8.23 |
| GA0 (hopper) | 11.59 / 11.65 | 8.61-9.61 | 20.09 / 21.54 | 18.65-24.77 |

kPa. Discharge values are means of 16 writes from 0.25 to 1.0 s after opening.

- At rest both models match the experiment at GA1 to GA5 and GA7. GA6 reads 14 % below
  it, just outside the 15 % band. The hopper wall reads 21 % high (11.6 against 9.6).
- In discharge the switch peak at GA0 falls inside the experiment's range for both
  models (20.1 and 21.5 kPa; P1's FEM 27.5). GA1 is within the band (-3 % MC, -7 % DP).
  GA2 is within it for MC (-12 %) but not for DP (-29 %). The experiment's discharge
  overpressure in the lower cylinder (GA1 to GA3, 1.4 to 1.9 times static) is only
  partly reproduced. Couto attribute it to density changes, which our model does not have.
- 12.5 mm (MC, 4 cores, 32 min) gives the same picture: GA0 19.3 kPa, GA1 5.6, GA2 4.2.
- The discharge rate is low: 46 kg/s (MC) and 40 kg/s (DP) at 25 mm, 69 kg/s at
  12.5 mm, against about 110 kg/s from Beverloo for a 0.35 m outlet. It rises with
  mesh refinement. Couto's real rate was throttled (0.6 to 4.9 kg/s), so it does not
  enter case A. It is the subject of H.4 and H.5.

## Decision

H.3 is closed with its done-when partly met: both models, at rest and in discharge,
match the experiment at most cells. The misses are recorded here and in STATUS. The
manual test stays strict and fails on GA6 and GA0 at rest (both models) and on GA2 in
discharge (DP). Package H.3b (DP robustness at partly filled and cut cells) is added.

## Consequences

- Validation runs use the T0 rule set with the wheat mu(I) values above, for both models.
- DP is not yet as robust as MC on free-surface discharge. MC is the recommended model.
