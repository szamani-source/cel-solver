# WP 3.7b: Silo fed from the top while it discharges

State: done
Branch: `claude/project-thread-53vkqx`

## Goal

The T0 tank discharges through its floor outlet while a solid inlet in its lid feeds
it (the inlet of WP 3.7a, decision 0060), as a test and as the first tutorial with
`inlets`.

## Read first

- `docs/decisions/0060-solid-inlet.md`
- `tests/outlet/outlet_checks.py` (the tank and outlet this reuses)
- `tutorials/siloDischarge` (the tutorial this copies)

## Build

- `tests/test_fed_silo.sh`, `tests/fedsilo/fed_silo_checks.py`.
- `tutorials/siloFed` (Allrun, plot.py, README).

## Done when

- `tests/test_fed_silo.sh` passes: mass balance with feed and discharge to 1e-10;
  discharge rate with feed R/2 and R within 10 % of the unfed R; matched feed keeps
  the tank mass within 3 %; no stress in empty cells.

## Out of scope

- Feeding on a wedge (WP 3.8c can reuse it).

## Log

- 2026-11-18: test written first and passing (about 2 min, E = 0.1 MPa, T0 rule
  set): unfed 73.7 kg/s/m; with feed 36.8 and 73.7 the rate out is -1.0 % and -2.7 %;
  matched feed holds the tank mass within 1.05 %; balance 4e-14. A first check 3
  (tank rate = out - feed) was dropped: it restates the mass balance. Tutorial
  `siloFed` run at 10 MPa for its plots.
