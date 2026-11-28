# celFoam tutorials

The tutorials are grouped in three categories, one folder each (decision 0140):

| Folder | Category | What is solved |
| --- | --- | --- |
| [`solid/`](solid/README.md) | Solid only | The granular or solid material and void; walls, outlets, inlets, bodies |
| [`solidFluid/`](solidFluid/README.md) | Solid and fluid | The above plus the gas in the pores and the void (version 2, decisions 0125-0135) |
| `solidFluidThermal/` | Solid, fluid and temperature | The above plus the energy equation (decision 0136) |

Every tutorial folder has `Allrun`, `Allclean` and a `README.md` with the same sections:
Setup, Numerics, How to run and Expected results, and the result figures it shows.
`tests/test_solid_tutorials.sh` checks the solid ones and `tests/test_gas_tutorials.sh`
the solid-fluid ones (layout, READMEs, a short run of each).

`hotSiloDischarge` (solid, gas and temperature) still sits at the top level until the
solid-fluid-temperature tutorials are gathered in their folder.

Copy the tutorials before running them: `celfoam tutorials ~/celfoam-cases` (installed
copy), or run them in place in a source checkout and `./Allclean` afterwards.
