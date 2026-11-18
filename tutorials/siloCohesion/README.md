# siloCohesion: the T0 tank with cohesion

Sia's target tank (T0, `tutorials/siloDischarge`) with a cohesive material: c = 500 Pa
instead of 0. The material settles for 0.3 s with the outlet closed, then the outlet
opens. With a narrow outlet the material arches over it and the flow stops; with a wide
one it flows. Allrun runs two outlet widths side by side.

| | |
| --- | --- |
| Tank | 0.40 m wide, filled to 0.60 m, outlet centred in the floor: 100 mm and 140 mm |
| Material | rho 1800 kg/m3, E 10 MPa (T0), nu 0.3, Mohr-Coulomb phi 35 deg, c 500 Pa (placeholder) |
| Tensile strength | c cot(phi) = 714 Pa (default; `sigmaT` sets another value) |
| Friction | mu(I) placeholders as `siloDischarge` (decision 0037) |
| Rule set | T0 rule set (decision 0050) |
| Walls, floor | friction 20 deg; outlet closed until 0.3 s; walls have no adhesion |
| Mesh | 10 mm square cells, 40 x 62, one cell deep |

Jenike's critical width for a slot outlet is B_c = sigma_c / (rho g), with the
unconfined strength sigma_c = 2 c cos(phi) / (1 - sin(phi)) = 1920 Pa here: B_c = 109 mm.
Below it a cohesive arch can carry the material above the outlet; above it, it cannot.

## How to run

```bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc   # tools/setup_openfoam.sh installs it
./Allwmake                                          # from the repository root, once
cd tutorials/siloCohesion
./Allrun                # outlets 100 and 140 mm, side by side (folders D100, D140)
./Allrun "60 100 140"   # other widths: multiples of 20 mm
python3 plot.py         # snapshots.png, mass.png, result.txt (needs numpy, matplotlib)
./Allclean
```

Change `c` (and optionally `sigmaT`) in `constant/materialProperties` before `./Allrun`
to try another cohesion; `E` 1e6 runs about three times faster and gives the same
critical width (decision 0054).

## What it shows

![Snapshots](snapshots.png)

The outlet opens at 0.3 s. With 100 mm (below Jenike's 109 mm, top row) a small piece
drops out of the opening and the rest stays put: the cohesive material arches over the
outlet. With 140 mm (bottom row) a flowing zone forms above the outlet and the tank
empties; grey cells inside the stream are void, where the material tears apart at its
tensile strength.

![Mass in the tank](mass.png)

| Outlet | Discharged in 1 s | Rate over the last 0.2 s |
| --- | --- | --- |
| 100 mm | 1.0 % of the fill | 4.3 kg/s/m (2.4 % of Beverloo with C = 1): stopped |
| 140 mm | 34 % of the fill | 171 kg/s/m (58 % of Beverloo with C = 1): flows |

The stopped tank still creeps slowly (a known issue, decision 0054). The test
`tests/test_cohesion.sh` brackets the critical width at c = 500 Pa between 100 and
120 mm (Jenike 109 mm), at E = 1 and 10 MPa; at c = 250 Pa it lies between 40 and 60 mm
(Jenike 55 mm).

## Dictionaries first used here

`c` was already an entry of the Mohr-Coulomb model; since WP 3.5 (decision 0054) it also
sets the tensile strength of the free-surface rules: contact faces between material
cells carry tension up to `sigmaT` (default c cot(phi)) and tear beyond it, and a cell
stretched beyond it opens void. With c = 0 (or `sigmaT 0`) the material separates as
soon as it pulls, as in `siloDischarge`.
