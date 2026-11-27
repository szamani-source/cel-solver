# stlSilo: the T0 tank built from STL walls

The siloDischarge tank (Sia's target problem T0) with its side walls and floor given as
STL surfaces instead of mesh boundaries (WP 4.3c, decisions 0047 and 0051). The walls
cut the 10 mm mesh at a quarter and three quarters of a cell, there is void above the
material, and the outlet is the mesh patch under the gap in the STL floor: closed while
the material settles for 0.3 s, then open for 2 s.

## Setup

| | |
| --- | --- |
| Tank | STL bodies `leftWall`, `rightWall`, `leftFloor`, `rightFloor` (constant/triSurface): 0.40 m wide (x 0 to 0.40), floor at y = 0.0125, gap 0.06 m (x 0.17 to 0.23) |
| Fill | 0.6 m above the floor (setFields), void above |
| Material | as siloDischarge: rho 1800, E 10 MPa, nu 0.3, phi 35 deg, c 0, mu(I) placeholders |
| Friction | 20 deg on every STL body (`friction` in constant/rigidBodies) |
| Outlet | mesh patch `outlet` under the gap, friction 20 deg until it opens at 0.3 s |
| Mesh | 46 x 65 cells of 10 mm, one cell deep; the STL tank does not follow it |

## Numerics

- Explicit Lagrange step plus remap onto the fixed mesh (decision 0013): MUSCL-Hancock
  face states with an acoustic Riemann solver, Jaumann stress rate (Hughes-Winget),
  time step from the CFL limit with `maxCo 0.4` (`system/controlDict`).
- Remap: unsplit, face based; default limiters (MC for the main fields, superbee for the
  plastic history, decision 0017).
- Material plus void with the free-surface rules (decisions 0021-0028); the T0 rule set:
  mu(I) friction with `shearCap`, `spreadLayers`, `flowImpedance` and `normalCap`
  (decision 0050).
- Settling: gravity ramps in over 0.1 s (`gravityRamp`) and the velocity is damped at
  250 1/s (`damping`) until the outlet opens (decision 0020). The outlet is a wall with
  friction until `openTime`, then an outflow face (decisions 0026, 0027).
- STL walls are immersed cut cells (decisions 0047, 0051): open areas and volumes per
  cell, the wall face carries the body's friction, cells cut to a small open volume
  are merged with a neighbour in the Lagrange step, covered cells hold ghosts.

## How to run

```bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc   # tools/setup_openfoam.sh installs it
./Allwmake                                          # from the repository root, once
cd tutorials/solid/stlSilo
./Allrun          # E = 10 MPa, about 9 minutes on one core
./Allrun 0.1      # softer material, about a minute
python3 plot.py   # snapshots.png, mass.png, rate.txt (needs numpy, matplotlib)
./Allclean
```

`alpha` is the material fraction of the open part of each cell (the part not inside a
body). Field `alphaBody` (written with the others) shows the walls in ParaView.

## Expected results

![Snapshots](snapshots.png)

![Mass in the tank](mass.png)

- No material leaves while the outlet is closed: the mass stays at 431.5 kg/m to 1e-15
  through the settling, so nothing leaks through the STL walls or floor.
- The flow looks like the mesh-wall tank: a fast stream over the gap, a slower
  converging zone, the top coming down nearly level.
- Steady rate 74.6 kg/s per metre of depth (Beverloo C = 0.90, run of 2026-11-27);
  the mesh-wall tank (`siloDischarge`) gave 60.2 in the same session. Single 10 mm
  rates scatter by about 15 % from round-off (decision 0046); in bands of three runs
  the STL tank reads 4.4 % above the mesh-wall tank (68.8 against 65.8 at 0.1 MPa,
  `tests/test_stl_silo.sh`, decision 0052).
