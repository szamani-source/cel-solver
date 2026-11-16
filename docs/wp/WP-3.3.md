# WP 3.3: C++ outlet boundary that opens

State: done (2026-11-15)
Branch: `claude/project-thread-v85uto`

## Goal

Add the outlet to the C++ tank: a closed frictional wall while the material settles,
then an outflow face (prototype rules 10 and 10b, decisions 0026 and 0027), with the
shear cap and rule 9 on (decision 0032). Check against the prototype at 10 mm.

## Read first

- `docs/decisions/0026`, `0027`, `0031`, `0032`; `prototype/void2d.py` (rule 10, 10b),
  `prototype/discharge2d.py`
- `src/celEulerian/eulerianSolid/eulerianSolidLagrange.C` (boundary faces),
  `eulerianSolidRemap.C` (`remapFreeSurface`)

## Build

- `materialProperties` `outlets { <patch> { friction <deg>; openTime <s>; } }`,
  `openOutlets()`, `discharged()`; rule 10b at open outlet faces.
- `celVoidCheck` writes the discharged mass in `history.dat`.
- `tests/test_outlet.sh`, `tests/outlet/outlet_checks.py`.

## Done when

`tests/test_outlet.sh` passes: nothing leaves before the opening time, material leaves
after it, tank plus discharged mass constant to 1e-10, no stress in void, rate C at least
0.8 and within 20 percent of the prototype at the same settings.

## Log

- 2026-11-15: Done in one session (decision 0033). Rate 74.1 kg/s/m against the
  prototype's 81.35 (-9 %); without rule 10b 27.4.
