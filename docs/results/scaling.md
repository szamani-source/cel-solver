# Parallel scaling (WP 6.3a)

Measured 2026-11-19 in the cloud container (4 cores), `tests/manual/test_scaling.sh`.
Case: `tutorials/siloDischarge` (T0 tank 0.4 x 0.6 m, E 10 MPa, T0 rule set) while it
settles (gravity ramp and damping, outlet closed): 0.05 s at 10 mm cells, 0.02 s at
5 mm; decomposition `simple`, split along x only.

| Cells | Processors | Cells | Steps | Wall time (s) | Cost (us per cell per step, x processors) | Speed-up | Efficiency |
|---|---|---|---|---|---|---|---|
| 10 mm | 1 | 2,480 | 1,082 | 10.9 | 4.05 | 1.00 | 100 % |
| 10 mm | 2 | 2,480 | 1,082 | 6.6 | 4.92 | 1.65 | 82 % |
| 10 mm | 4 | 2,480 | 1,082 | 4.1 | 6.19 | 2.62 | 65 % |
| 5 mm | 1 | 9,920 | 866 | 35.9 | 4.18 | 1.00 | 100 % |
| 5 mm | 2 | 9,920 | 866 | 19.5 | 4.53 | 1.84 | 92 % |
| 5 mm | 4 | 9,920 | 866 | 10.9 | 5.07 | 3.30 | 82 % |

Every run took the same number of steps; the mass agrees to 1e-13 (sums in a different
order).

## Reading it

- A cell step costs about 4 microseconds on one core with the full T0 rule set (mu(I),
  flow impedance, caps, void), against 2.2 for the bare elastic cycle (WP 2.4b).
- Efficiency falls when a processor holds few cells: 620 cells each (10 mm on 4) gives
  65 %, 2,480 each (5 mm on 4) 82 %. Keep at least a few thousand cells per processor.
- Mesh refinement costs 8 times per simulated second in 2D (4 times the cells, half the
  step), 16 times in 3D. A 5 mm T0 tank discharge (2.3 s) takes about 70 min on one
  core and about 20 min on four.
- Not measured: more than 4 cores, 3D meshes, bodies (the hex8 subcycles run on every
  processor in full, so they do not speed up), or scotch decomposition (simple along x
  splits the empty void above the fill evenly with the material, which is fine here).
