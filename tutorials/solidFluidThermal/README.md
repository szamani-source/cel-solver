# Solid, fluid and temperature tutorials

The material, the gas in its pores and the energy equation (version 2, decision 0136).
Each folder has `Allrun`, `Allclean` and a `README.md` with the sections Setup,
Numerics, How to run and Expected results (decision 0140).

| Tutorial | What it shows | Run time |
| --- | --- | --- |
| `hotSiloDischarge` | Hot sand under cold sand in the T0 tank: the temperature carried by the discharge, cooled side walls, the outlet temperature over time | about 10 min (1 core), 1 min with `./Allrun 0.1` |

`tests/test_thermal_tutorials.sh` checks the layout and the READMEs;
`tests/test_hot_silo_tutorial.sh` runs `hotSiloDischarge` to 0.8 s.
