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
- 2026-11-22: M.1 to M.9 done in the same session and the same PR as M.0 (CI takes about
  3 h). The FAQ and known issues page (planned for M.10) was written with them. 43
  references; doi.org is blocked, so the DOIs of the less familiar entries were checked
  by web search, which corrected Pailha and Pouliquen (2009) to 10.1017/S0022112009007460.
  Miller and Colella (2001) has no DOI link (not confirmed); Jop et al. (2005) links to
  arXiv. Run `tools/check_manual.py --online` where doi.org is reachable. Fixed for phone
  width: KaTeX's hidden MathML widened the page (`.katex { position: relative }`), and
  long inline paths now wrap. Next: M.10 preview and review.
