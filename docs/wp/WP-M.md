# WP M.0 to M.10: HTML user manual

State: in progress
Branch: `claude/project-thread-lsis47` (thread "Build HTML user manual for solver")

## Goal

A comprehensive user manual in `docs/manual/`: install, quick start, case setup
reference, meshing guide, theory with citations, validation, references, FAQ.
Plan and ETAs approved by Sia 2026-11-22 (about 25 h). Decision 0116.

## Read first

- `docs/decisions/0116-html-user-manual.md` (format and checks)
- Sources: `/mnt/project-files/notes/celFoam-equations.md` (theory, as of 2026-11-16),
  `docs/constitutive-models.md`, `docs/decisions/`, `docs/results/validation-report.md`,
  `docs/spec/validation-references.md`, `tutorials/*/README.md`

## Build

- One page per package; add it to `NAV` in `docs/manual/assets/manual.js` and its
  required terms to `REQUIRED` in `tools/check_manual.py`.

## Done when

- `tests/test_manual.sh` passes with the package's page present (see the M table in
  `docs/work-packages.md`).

## Out of scope

- New simulations (Sia: "no runs"). Figures reuse existing results.

## Log

- 2026-11-22: M.0 done. Skeleton (CSS, JS menu, vendored KaTeX), home page, reference
  list, `tools/check_manual.py` and `tests/test_manual.sh`. doi.org and Crossref are
  blocked from the cloud container. Next: M.1 installation.
