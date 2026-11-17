# WP 3.6s: mesh dependence of the collapse shear band (split from 3.6b)

State: done (finding recorded; no fix in this package)
Branch: `claude/project-thread-in8o8d`

## Goal

Find out whether the shear bands of the C++ column collapse depend on the mesh, with the
plain T0 Mohr-Coulomb and with mu(I) (decision 0037 says T0's Jop law is not
regularised). Split from WP 3.6b on 2026-11-17 so the inclined layer stays its own package.

## Read first

- `docs/decisions/0032-*.md` (C++ collapse), `docs/decisions/0037-*.md` (mu(I)).
- `tests/collapse/collapse_checks.py` (set-up reused).

## Build

- `tests/collapse/shear_band.py`, `tests/manual/test_shear_band.sh` (about 9 min).
- Overflow guard in `FrictionLaw::inertialNumber` (C++) and `inertial_number` (Python).

## Done when

- `tests/manual/test_shear_band.sh` runs a = 0.5, E = 1 MPa on 10, 5 and 2.5 mm cells with
  both laws, all runs finish with mass kept to 1e-10, and the band widths and run-outs are
  recorded here and in decision 0038. Met 2026-11-17.

## Out of scope

- The cure (regularisation or a length scale): a later package, decision 0038.
- The inclined layer (WP 3.6b).

## Results (2026-11-17)

Collapse a = 0.5 (L_i 0.2 m, H_i 0.1 m), E = 1 MPa, shear cap and rule 9 on, floor friction
35 deg, t = 1 s. Band width = horizontal full width at half the peak shear rate, rows at 20
to 60 % of H_i, averaged. Run-out with the front where the layer is 5 mm thick (in brackets
the test_collapse.sh measure, half a cell of mass, which itself depends on the cell).

| Law | Cell | Band width 0.1 s / 0.2 s | Peak shear rate 0.1 / 0.2 s | Run-out 0.5 s / 1 s | Height at x = 0.2 m, 1 s | Wall time |
| --- | --- | --- | --- | --- | --- | --- |
| MC | 10 mm | 25 / 45 mm (2.5 / 4.5 cells) | 10.3 / 2.9 1/s | 0.20 / 0.30 (0.20 / 0.30) | 62 mm | 10 s |
| MC | 5 mm | 5.6 / 8.8 mm (1.1 / 1.8 cells) | 22.7 / 15.5 1/s | 0.25 / 0.33 (0.28 / 0.33) | 56 mm | 65 s |
| MC | 2.5 mm | 8.1 / 7.5 mm (3.2 / 3.0 cells) | 29.7 / 22.6 1/s | 0.40 / 0.43 (0.56 / 0.58) | 50 mm | 466 s |
| mu(I) | 10 mm | 18 / 38 mm (1.8 / 3.8 cells) | 10.4 / 3.2 1/s | 0.20 / 0.25 (0.20 / 0.30) | 62 mm | 11 s |
| mu(I) | 5 mm | 5.6 / 12.5 mm (1.1 / 2.5 cells) | 22.8 / 11.6 1/s | 0.25 / 0.33 (0.28 / 0.33) | 56 mm | 67 s |
| mu(I) | 2.5 mm | 8.3 / 9.5 mm (3.3 / 3.8 cells) | 30.1 / 19.1 1/s | 0.44 / 0.56 (1.01 / 0.75) | 47 mm | 521 s |

Mass kept to 1e-13 in every run. Top height 90 to 92 mm on every mesh.

What it says:

- **The band is mesh dependent.** It is 1 to 4 cells wide on every mesh, and the peak shear
  rate rises with refinement (10, 23, 30 1/s at 0.1 s). Width in mm does not settle
  (25, 6, 8 mm at 0.1 s). This is the usual localisation of a rate-independent softening-free
  frictional model with non-associated flow (psi = 0), so the band width is set by the cell.
- **mu(I) does not cure it.** At these rates I is about 0.05 (mu about 0.04 above mu_s),
  so mu(I) and plain MC give the same bands at 10 and 5 mm; at 2.5 mm mu(I) differs only in
  the late run-out (a thin layer keeps sliding at 1 s). As decision 0037 feared, the plain
  Jop law adds no length scale.
- **Run-out does not converge** (5 mm front: 0.30, 0.33, 0.43 for MC at 1 s); the height
  behind the front does (62, 56, 50 mm at x = 0.2 m, steps shrinking roughly by half), and
  the top does. The 10 mm run-outs in `tests/test_collapse.sh` are therefore mesh-bound
  numbers, not converged ones.
- **Bug found and fixed:** at 5 mm the mu(I) run stopped with a floating-point trap in
  `inertialNumber`: a free-surface cell with mean pressure about 1e-320 Pa made rho_s/P
  overflow. Such P now counts as P = 0 (I infinite, mu = mu_2), in C++ and Python alike;
  `tests/test_mu_i.sh` still passes.

## Log

- 2026-11-17: study run and recorded (table above), decision 0038. Next: WP 3.6b (inclined
  layer, its own thread); the band cure is listed in decision 0038 and STATUS as a later
  package.
