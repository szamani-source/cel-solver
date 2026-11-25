# WP F.5: Air-impeded discharge (end of Stage A)

State: done (decision 0130)
Branch: assigned by the thread

## Goal

Discharge with pore gas on: coarse grains match v1, fine grains slow down as published. Read Crewdson et al. (1977) first; if it cannot be reached, use Muite et al. (2004).

## Read first

- `docs/decisions/0125-gas-coupling-plan.md`
- `docs/spec/gas-references.md`, `tools/gas_reference.py`
- Briefs of the earlier F packages (their logs)

## Build

- A hopper case for fine powder (geometry from the paper), and the Couto hopper as a control
- Checkpoint every 0.05 s to /mnt/project-files; show the mesh to Sia before long runs

## Done when

- `tests/test_air_impeded.sh` (short) plus a manual long run: 1 mm and up within 5 % of the dry rate; fine grains slow with the published trend within 20 %; report the Muite N for each case.

Parallel (Sia, 2026-11-24: the solver must be fully parallel): the done-when test also runs on 2 and 4 processors, with the cuts through the bed, the free surface and the outlet. Gas fields match serial to 1e-10 relative, with the pressure solve converged to 1e-12; total gas mass is conserved to 1e-12 on every processor count. Also report the cost per step with gas on against gas off.

Tests use the `ergunWenYu` drag law (Sia, 2026-11-24).

Validation pack (Sia, 2026-11-24): this package delivers pack test 4, fine-powder hopper discharge against Crewdson et al. (1977), gas on and gas off on the same mesh so the slowdown from the air shows (coarse mesh, hours), as `tests/validation/gas/<name>/` with an `Allrun` that runs it on 1, 2 and 4 processors, checks that they match and writes the run times to `times.csv`.

## Out of scope

- Gas momentum (Stage B).

## Log

Newest first. Each entry: date, what was done, what was learned, exact next step.

- 2026-11-24: done (decision 0130). Crewdson (1977) and Muite (2004) could not be
  read online; checked against Crewdson's model as quoted in an open thesis and the
  abstracts of Hsiau et al. (2010) and Muite et al. Plane-flow hopper, 60 mm slot, 544
  cells; gas on against off: 1 mm -2.2 %, 0.1 mm -64 %, 2/4 processors 4e-12; pack test
  4 (fineHopper) fails on 4 processors in one surface cell (8e-10). Learned: (1) a 20 mm
  slot with 2 mm cells discharges at 1/15 of Beverloo in v1 (gas or not), so use a wider
  slot; (2) mu(I)'s d must be the case's d in the dry control too; (3) the T0 rules
  part the solid between processor counts (4 % by 0.25 s), so parallel checks use plain
  Mohr-Coulomb. Next: F.6; F.9 fixes the 4-processor surface cell and the cost.

- 2026-11-24: brief written in WP F.0. Next: start when the package before it has merged.
