# WP 3.7a: Solid inlet that feeds material into the domain

State: done
Branch: `claude/project-thread-53vkqx`

## Goal

A patch can feed fresh solid material into the domain (something Abaqus CEL cannot
do, Sia 2026-11-18): set speed over time, density, stress, no plastic history, mass
counted.

## Read first

- `docs/decisions/0033-outlet-cpp.md` (outlets, the pattern this follows)
- `docs/decisions/0031-free-surface-cpp.md` (slab rules)
- `src/celEulerian/eulerianSolid/eulerianSolidRemap.C` (`remapFreeSurface`)

## Build

- `inlets` in materialProperties, read in `eulerianSolid.C` (`readInlets`,
  `updateInlets`, `charged()`).
- Inlet face in the Lagrange step boundary loop; inflow volume, mass and inlet values
  in `remapFreeSurface` and `remapMassWeighted`.
- `celVoidCheck` history gets the fed mass (column 23); celFoam prints discharged and
  fed mass.

## Done when

- `tests/test_inlet.sh` passes: closed inlet lets nothing in; fed mass exact (1e-12);
  box mass = fed mass (1e-12); no stress in empty cells; fed density kept; the deposit
  comes to rest; 2 processors match serial.

## Out of scope

- A tutorial with an inlet and the fed silo (WP 3.7b).
- Inlets without `freeSurface`, inlets on STL bodies, outflow through an inlet.

## Log

- 2026-11-18: done. Test written first, then the inlet (decision 0057). Counting an
  open inlet as full for the slab side made the stream fall in place; as empty it
  falls and lands. Plain MC packs the landing (+4.3 %); the test uses the T0 rule set
  (+0.82 %). All seven checks pass in 25 s; test_void and test_outlet unchanged. Next:
  WP 3.7b.
