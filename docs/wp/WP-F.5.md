# WP F.5: Air-impeded discharge (end of Stage A)

State: not started
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

Tests use the `ergunWenYu` drag law (Sia, 2026-11-24).

## Out of scope

- Gas momentum (Stage B).

## Log

Newest first. Each entry: date, what was done, what was learned, exact next step.

- 2026-11-24: brief written in WP F.0. Next: start when the package before it has merged.
