# WP 3.6d: silo discharge with mu(I) and the flow impedance

State: closed, done-when not met (2026-11-17)
Branch: `claude/project-thread-8hqs9j`

## Goal

Run Sia's target T0 (tank, settle, open, discharge) with mu(I) and `flowImpedance yes`
at E = 0.1, 1, 10 MPa; check the rate is independent of E; compare with Beverloo and the
silo references; package it as the first case for Sia to look at.

## Read first

- Decisions 0034, 0037, 0040, 0044, 0045; `tests/discharge/discharge_checks.py`
- `tutorials/siloDischarge/README.md`

## Done when

`tests/manual/test_discharge_stiffness.sh` (now with mu(I)) passes at 10 mm and 5 mm, C
within decision 0019.

## Log

- 2026-11-17: One session. First 10 mm runs: 62 / 50 / 12 kg/s/m. Crash at step 1
  (P ~ 1e-320 in flowViscosity) fixed. Bisected the flow impedance with temporary
  switches: the nearly-full snap of decision 0044 inflated the stream above the outlet
  at 10 MPa (density 967). Bounded the snap by the stress-implied density: 65.5 / 53.8 /
  67.7 at 10 mm, 77.6 / 61.2 at 5 mm (0.1 / 1 MPa); layer still passes. 1 MPa stays
  about 20 % slow; 10 mm rates scatter about 15 % from round-off. Made the case a celFoam
  tutorial (`tutorials/siloDischarge`, plots). Not run: 5 mm at 10 MPa (about 70 min).
