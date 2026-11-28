# Solid-only tutorials

Cases with the solid (granular) material, void, walls, outlets, inlets and rigid
bodies: no gas and no temperature. Each starts from a validated or regression-tested
case, so its numbers can be trusted as far as that test goes (column "Checked by").
Start with `smokeBox` (seconds), then `siloDischarge`.

| Tutorial | What it shows | Run time | Checked by |
| --- | --- | --- | --- |
| `smokeBox` | 10 x 10 x 10 elastic box with a velocity jump: checks the install | seconds | `tests/test_smokeBox.sh` |
| `siloDischarge` | Sia's target tank (T0): settle, open the floor outlet, discharge; Mohr-Coulomb with mu(I) | 11 min (`./Allrun 0.1`: 40 s) | `tests/test_discharge.sh`, Beverloo |
| `siloCohesion` | The T0 tank with cohesion: arching over a narrow outlet, flow through a wide one | 7 min | `tests/test_cohesion.sh`, Jenike critical width |
| `siloFed` | The T0 tank fed through a solid inlet in the lid while it discharges | 30-60 min (`./Allrun 0.1`: minutes) | `tests/test_fed_silo.sh` |
| `stlSilo` | The T0 tank with walls and floor from STL files (immersed cut cells) | 10 min | `tests/test_stl_silo.sh` |
| `stlHopper` | A 60 deg hopper from STL walls, mass flow to empty | 2 min | `tests/test_stl_silo.sh` |
| `coutoHopper` | The validated Couto wheat silo on the axisymmetric wedge, wall-aligned blockMesh mesh with an STL shell; Mohr-Coulomb or Drucker-Prager | 2 h (standard) or 8 h (fine) on 4 cores | Couto et al. 2013, Zheng and Yu 2015, Beverloo |
| `setupPlaneTank` | The same kind of tank written by `celfoam setup` from a short YAML file | about 1 h (`./Allrun 0.1`: minutes) | `tests/test_celsetup.sh` |
| `stripFooting` | Coupled Eulerian-Lagrangian: a rigid STL footing pushed into undrained clay | 1 min (`--fine`: 9 min) | `tests/test_footing.sh`, Prandtl 2 + pi |
| `rigidBodyMotion` | Prescribed rigid-body kinematics: a fixed wall, a pushed and a spinning cube | seconds | `tests/test_rigid_body.sh` |

Features by tutorial: mesh walls and outlets (`siloDischarge`), cohesion and tensile
strength (`siloCohesion`), solid inlets (`siloFed`), STL walls (`stlSilo`, `stlHopper`,
`coutoHopper`), axisymmetric wedge and angled hex cells (`coutoHopper`), Drucker-Prager
(`coutoHopper dp`), moving bodies and the force on them (`stripFooting`,
`rigidBodyMotion`), the case generator (`setupPlaneTank`).

Not a tutorial yet: the compaction cap (decision 0113) is checked at one point and in a
box (`tests/test_cap.sh`, `tests/test_cap_stl.sh`) but did not close the at-rest gap in
the Couto silo, so no silo tutorial uses it.

Every README has the sections Setup, Numerics, How to run and Expected results. The
figures beside each README come from full runs of the tutorial as shipped (dates in
each README); `tests/test_solid_tutorials.sh` runs a short stretch of every case.
