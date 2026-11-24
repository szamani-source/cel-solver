# Gas validation pack test 3: pore-pressure decay against Terzaghi

Fine quartz powder d 0.1 mm, gas fraction 0.4, E 1 MPa, nu 0.3 (M 1.346 MPa); column 1.0 m (1 x 80 cells), drained top loaded with 2000 Pa at t = 0; air at 293.15 K; drag ergunWenYu; damping 50 1/s. c = 0.1393 m^2/s, p_i = 316.75 Pa. Written by `terzaghi_column.py`.

Serial run against the series:

| T | U | Terzaghi | error | profile, largest error / p_i |
| --- | --- | --- | --- | --- |
| 0.05 | 0.2511 | 0.2523 | -0.50 % | 0.008 |
| 0.1 | 0.3545 | 0.3568 | -0.66 % | 0.006 |
| 0.2 | 0.5012 | 0.5041 | -0.57 % | 0.005 |
| 0.5 | 0.7614 | 0.7640 | -0.34 % | 0.004 |
| 1 | 0.9299 | 0.9313 | -0.15 % | 0.002 |

| processors | steps | wall [s] | worst U error | vs serial | gas mass | solves per step |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 41741 | 20.9 | 0.66 % | 0e+00 | 1e-15 | 2.00 |
| 2 | 41741 | 27.1 | 0.66 % | 1e-12 | 9e-16 | 2.04 |
| 4 | 41741 | 34.0 | 0.66 % | 9e-13 | 5e-16 | 2.03 |

All checks pass: U within 3 % of Terzaghi, profiles within 0.03 p_i from T = 0.1, processor counts agree to 1e-10 of p_i, gas mass to 1e-12. Before T = 0.05 the elastic waves of the sudden load (not in the quasi-static theory) are still being damped.
