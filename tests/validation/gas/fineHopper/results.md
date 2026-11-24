# Gas validation pack test 4: fine-powder hopper, gas on and off

Plane-flow hopper, slot 60 mm, walls 30 deg from vertical, 544 body-fitted cells; quartz powder d 0.1 mm, gas fraction 0.4, plain Mohr-Coulomb phi 30, psi 1, E 1 MPa; drag ergunWenYu; outlet opens at 0.15 s, rate from 0.3 to 0.5 s. Written by `fine_hopper.py`.

Gas on against gas off (serial): **16.92 against 32.94 kg/(m s), -48.6 %** (limit: at least 20 % slower). Crewdson's sqrt(1 + G/(rho_b g)) from the run's outlet gradient: 0.390. Muite's N = U_t/U_out: 3.01 (the gas matters below about 10).

| run | processors | steps | wall [s] | rate [kg/(m s)] | vs serial | gas mass | solves per step |
| --- | --- | --- | --- | --- | --- | --- | --- |
| dry | 1 | 12229 | 15.4 | 32.937 | 0e+00 | 0e+00 | 0.00 |
| dry | 2 | 12229 | 10.6 | 32.937 | 3e-12 | 0e+00 | 0.00 |
| dry | 4 | 12229 | 7.8 | 32.937 | 3e-12 | 0e+00 | 0.00 |
| gas | 1 | 12208 | 36.4 | 16.921 | 0e+00 | 4e-15 | 8.77 |
| gas | 2 | 12208 | 30.2 | 16.921 | 7e-11 | 2e-15 | 8.77 |
| gas | 4 | 12208 | 25.1 | 16.921 | 8e-10 | 2e-15 | 8.77 |

Cost per step, gas on against off: 1 processor(s) 2.36x, 2 processor(s) 2.85x, 4 processor(s) 3.24x (F.9 target 1.5x).

- gas on 4 processors: U differs from serial by up to 8e-10 in cells (index, alpha) [(449, 0.48), (449, 0.58), (449, 0.68), (449, 0.78), (449, 0.88), (449, 0.98)]; full cells within 6e-11

All checks FAIL.
