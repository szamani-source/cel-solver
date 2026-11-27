# WP TU.1: Solid-only tutorials

State: done
Branch: `claude/project-thread-flrj1e`

## Goal

A set of ready-to-run solid-only tutorials, built from the validated cases, each with a
README (setup, numerics, how to run, expected results) and result figures the HTML
manual can use (Sia, 2026-11-26: tutorials in three categories, then the manual).

## Read first

- `tutorials/` (all), `docs/decisions/0117`, `0107`, `0055`, `0054`, `0060`.

## Build

- `tutorials/README.md`, `tutorials/solid/README.md`, the solid tutorials moved into
  `tutorials/solid/`, new `stripFooting`, README rewrites, plot scripts, figures.
- Path updates in tests, `install.sh`, `tools/celfoam`, CI install workflow, INSTALL.md,
  README.md, `docs/manual/quickstart.html`.
- Decision 0140.

## Done when

- `tests/test_solid_tutorials.sh` passes: every solid tutorial listed in the index, with
  executable Allrun and Allclean, the four README sections and the figures it shows;
  short runs of every case end cleanly (about 1.5 min).

## Out of scope

- Gas and temperature tutorials (`solidFluid/`, `solidFluidThermal/`), the manual
  rebuild, solver changes.

## Log

- 2026-11-27: Layout, test, READMEs, full runs for the figures (on 4 cores side by
  side, so run times in the READMEs are single-run times from earlier work packages
  where they were measured alone). Couto figures from the 4 mm validation runs.
