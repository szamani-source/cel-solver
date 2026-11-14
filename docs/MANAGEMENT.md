# How this project is managed

Purpose: let Claude pick up any part of the work, weeks later and with no memory of
earlier conversations, and know within five minutes what exists, why, and what to do next.
Sia can use the same files to see progress without reading code.

## The rule everything follows

**If it is not written in the repository, it does not exist.** Conversations are
forgotten; files, tests and git history are not.

## The five pieces

| Piece | File | Answers | Updated |
| --- | --- | --- | --- |
| Map | `CLAUDE.md` | What is this code and how is it built? | When the architecture changes |
| Board | `docs/STATUS.md` | What is done, in progress, next, and blocked on Sia? | End of every session |
| Briefs | `docs/wp/WP-x.y.md` | What exactly does this package do, and what happened so far? | During the package |
| Decisions | `docs/decisions/NNNN-*.md` | Why did we choose this? | When a choice is made |
| Tests | `tests/test_*.sh` | Does finished work still work? | With every package |

## Work packages

- The whole project is a list of small packages in `docs/work-packages.md`.
- A package is sized for one to three sessions and has exactly one "done when" test.
- Each package has a brief (`docs/wp/WP-x.y.md`, from `docs/wp/TEMPLATE.md`), written
  before work starts. The brief lists what to read first, what to build, what is out of
  scope, and ends with a dated log.
- A package that grows too large is split, and the split is recorded.
- One package = one project thread = one git branch = one pull request.

## Session routine

Start (about five minutes):

1. Read `CLAUDE.md`.
2. Read `docs/STATUS.md` and pick the package marked "next" (or the one Sia names).
3. Read that package's brief, including its log, and only the files the brief lists.
4. Run `tests/run_tests.sh` to confirm the starting point is green.

End (always, even if the package is not finished):

1. Add a dated entry to the brief's log: what was done, what was learned, exact next step.
2. Update `docs/STATUS.md`.
3. Write a decision note for any choice made.
4. Commit and push. Open or update the pull request `WP x.y: <title>`.

## What Sia does

- Answers items under "Waiting on Sia" in `docs/STATUS.md`.
- Reviews and merges pull requests (or tells Claude to merge).
- Judges physics results at each phase gate.

## What keeps this honest

- A package is "done" only when its test passes in CI, never by assertion.
- Existing tests are never weakened to get a green run.
- Anything not verified is written as "not verified" in the log and the status board.
