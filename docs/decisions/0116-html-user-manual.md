# 0116: User manual as static HTML pages in the repository

Date: 2026-11-22. Status: accepted (Sia approved the plan, 2026-11-22 15:01 UTC).

## Context

Sia asked for a comprehensive user manual: installation, case setup, meshing, the
theory with citations to published sources, and the validation cases, as a series of
HTML pages. It must stay correct as the code changes and must work on Sia's laptop and
in the Linux package, where there may be no internet.

## Decision

- Hand-written static HTML in `docs/manual/`, one page per topic, no build step.
- One stylesheet (`assets/manual.css`) and one script (`assets/manual.js`). The script
  holds the page list (`NAV`) and builds the side menu, on-page contents and
  previous/next links, so adding a page means one line in `NAV`.
- Equations in LaTeX between `$$ ... $$` (display) and `\( ... \)` (inline), typeset by
  KaTeX 0.16.11 vendored in `assets/katex` (MIT licence, woff2 fonts only, 0.6 MB).
  Single `$` is not a delimiter, so shell commands are safe.
- Citations are links `<a class="cite" href="references.html#ref-<key>">`; one reference
  list for the whole manual.
- `tools/check_manual.py`, run by `tests/test_manual.sh` in CI, checks the page list,
  links, anchors, citations both ways, tag balance, equations (typeset in strict mode
  when node is installed) and per-page required terms. Once `case-setup.html` exists it
  also fails if the C++ code reads a dictionary keyword the page does not document, so
  new options cannot go undocumented.
- Figures come from existing results; the manual work runs no simulations (Sia).
- Packages M.0 to M.10 in `docs/work-packages.md`, one pull request each.

## Consequences

- The manual opens from disk (`docs/manual/index.html`) and on any web server.
- DOIs cannot be checked from the cloud container (doi.org and Crossref are blocked by
  the network policy); `check_manual.py --online` does it where the network allows, and
  M.9 checks them by web search otherwise.
- Every new keyword in the code needs a line on the case setup page, or CI fails.
