# Thermal validation pack

Written by `thermal_pack.py` (decision 0136).

| case | processors | wall [s] | result | error | limit | vs serial |
| --- | --- | --- | --- | --- | --- | --- |
| schumann | 1 | 12 | theta_g 0.38 %, theta_s 0.26 % | 0.38 % | 1 % | 0e+00 |
| schumann | 2 | 14 | theta_g 0.38 %, theta_s 0.26 % | 0.38 % | 1 % | 3e-13 |
| schumann | 4 | 14 | theta_g 0.38 %, theta_s 0.26 % | 0.38 % | 1 % | 4e-13 |
| graetz | 1 | 101 | Nu 7.4570 against 7.5407 | 1.11 % | 3 % | 0e+00 |
| graetz | 2 | 62 | Nu 7.4570 against 7.5407 | 1.11 % | 3 % | 3e-11 |
| graetz | 4 | 45 | Nu 7.4570 against 7.5407 | 1.11 % | 3 % | 3e-11 |
| movingBed | 1 | 45 | h 291.3 against 292.8 W/(m^2 K), local Nu 0.45 % | 0.52 % | 5 % | 0e+00 |
| movingBed | 2 | 28 | h 291.3 against 292.8 W/(m^2 K), local Nu 0.45 % | 0.52 % | 5 % | 2e-14 |
| movingBed | 4 | 18 | h 291.3 against 292.8 W/(m^2 K), local Nu 0.45 % | 0.52 % | 5 % | 2e-14 |

All checks pass.
