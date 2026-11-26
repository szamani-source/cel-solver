# Gas validation pack test 2: onset of fluidisation against Wen and Yu

Quartz sand d 0.2 mm, packed at gas fraction 0.38, bed 0.05 m deep (8 x 20 cells, 20 void rows above), air at 293.15 K ramped to 0.0482 m/s over 0.28 s; model gasFlow with the pore gas momentum, drag ergunWenYu. Drag law's own U_mf at this packing 0.03355 m/s. Written by `fluid_onset.py`.

| processors | steps | wall [s] | U_mf [m/s] | Wen-Yu [m/s] | error | lift from [m/s] | vs serial to the onset | after (report) | gas mass |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 25370 | 70.7 | 0.03341 | 0.03446 | -3.05 % | 0.03637 | 0e+00 | 0e+00 | 8e-16 |
| 2 | 25370 | 71.2 | 0.03341 | 0.03446 | -3.05 % | 0.03637 | 1e-11 | 2e+00 | 6e-16 |
| 4 | 25370 | 73.8 | 0.03341 | 0.03446 | -3.05 % | 0.03637 | 1e-11 | 1e+00 | 8e-16 |

All checks pass: within 10 % of Wen and Yu, processor counts agree to 1e-10 of each field's largest value (p, Ug, U) up to the onset, gas mass to 1e-12. After the onset the lifting bed parts serial from parallel (round-off amplified; reported, F.9, decision 0132).
