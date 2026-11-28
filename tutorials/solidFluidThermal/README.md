# Solid, fluid and temperature tutorials

The material, the gas in its pores and in the void, and the energy equation (version 2,
decision 0136). Each folder has `Allrun`, `Allclean`, `plot.py` and a `README.md` with
the sections Setup, Numerics, How to run and Expected results (decision 0140). The first
three are the validated cases of the temperature series, set up as templates.

| Tutorial | What it shows | Checked against | Run time (1 core) |
| --- | --- | --- | --- |
| `hotAirPackedBed` | Hot air blown through a cold packed bed: gas and grain temperatures, grain-gas exchange (`twoTemperature`, Gunn) | Schumann (1929): 0.5 % | 15 s |
| `movingBedHeatedWall` | Particles flowing in plug flow past a plate heated by a uniform flux (`equilibrium`, wall `q`) | Adapa et al. (2024): mean h -0.5 % | 45 s |
| `heatedAirChannel` | Laminar air cooled between two plates in the void (`gasFlow`, `twoTemperature`, fixed wall temperature) | Shah and London (1978), Nu 7.5407: -1.1 % | 100 s |
| `hotSiloDischarge` | Hot sand under cold sand in the T0 tank: temperature carried by the discharge, cooled side walls | energy balance, outlet temperature | 10 min, 1 min with `./Allrun 0.1` |

`tests/test_thermal_tutorials.sh` checks the layout and the READMEs and runs a short
stretch of each case; `tests/test_hot_silo_tutorial.sh` runs `hotSiloDischarge` to 0.8 s.
