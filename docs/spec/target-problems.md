# Target problems, units and material data (WP 0.3)

Status: AGREED for the primary problem (Sia, 2026-11-14). Numbers marked "placeholder"
were chosen by Claude and can be changed without changing the plan.

## Units

SI throughout: m, kg, s, Pa. Stress positive in tension. Gravity 9.81 m/s^2, acting in
-y for 2D cases (x-y plane, one cell thick in z) and in -z for 3D cases.

## T0 (primary): 2D discharge from a rectangular tank through a bottom opening

Sia's description (2026-11-14): a 2D rectangular tank with a small opening at the
bottom. The simulation starts with the outlet closed; once the material is initialised
(settled under gravity), the opening opens and the material exits. Start with zero
cohesion; cohesion is added later.

| Item | Value | Status |
| --- | --- | --- |
| Dimension | 2D plane strain | Agreed |
| Tank | Rectangle, width W = 0.40 m, walls rigid | Placeholder |
| Fill | Height H = 0.60 m of material, void above (tank height 0.80 m) | Placeholder |
| Outlet | Centred in the flat bottom, width D = 0.06 m (also 0.04 and 0.08 m for the rate law) | Placeholder |
| Stage 1 | Outlet closed; gravity applied; material settles to rest | Agreed |
| Stage 2 | Outlet opens at a set time; material leaves the domain through it | Agreed |
| Material | Cohesionless frictional (Mohr-Coulomb, c = 0), see below | Agreed in kind |
| Walls | Coulomb friction, wall friction angle 20 deg; smooth walls as a variant | Placeholder |
| Mesh | Uniform, 5 mm cells (outlet 12 cells wide), then 2.5 mm | Placeholder |
| Later | Cohesion c > 0 (arching, rat-holing); hopper (inclined) bottom; 3D | Planned |

This problem needs no Lagrangian body: the tank is the mesh boundary. It needs gravity,
frictional walls, a free surface (void above the material) and a boundary that switches
from wall to outflow. The plan is re-ordered around it (decision 0016).

What "correct" means for T0, without a commercial code (reference values are fixed in
WP 3.0, with sources):

- Settled state: wall and bottom pressures follow Janssen's saturation with depth, not
  the hydrostatic line.
- Discharge rate is nearly constant in time and independent of the fill height while
  the bed is deeper than about one tank width (unlike a liquid).
- Discharge rate scales as `rho sqrt(g) D^(3/2)` per unit depth (Beverloo/Hagen form
  in 2D) across the three outlet widths.
- Flow pattern: a central flowing channel above the outlet, stagnant zones in the bottom
  corners, a crater forming in the free surface.
- Mass is conserved: mass in the tank plus mass through the outlet stays constant.

## Other targets (later, unchanged)

| # | Problem | Needs |
| --- | --- | --- |
| T1 | Rigid strip footing on undrained clay (Prandtl, Nc = 5.14) | Rigid Lagrangian body, contact |
| T2 | Cone or pile penetration | Rigid Lagrangian body, contact |
| T3 | Projectile into soil | Deformable FE body |

## Material data

| Material | Model | Parameters | Status |
| --- | --- | --- | --- |
| Bulk solid for T0 | Mohr-Coulomb, c = 0, tension cut-off 0 | rho 1800 kg/m3, E 10 MPa, nu 0.3, phi 35 deg, psi 0 deg, c 0, sigma_t 0 | Agreed by Sia (2026-11-14): clay placeholder density and stiffness, friction angle 35 deg, zero cohesion to start |
| Same, cohesive (later) | Mohr-Coulomb | as above with c 1 to 5 kPa | Planned |
| Undrained clay (T1, T2) | Tresca (von Mises with su) | rho 1800 kg/m3, E 10 MPa, nu 0.49, su 20 kPa | Placeholder |
| Steel (T3) | von Mises, later Johnson-Cook | rho 7850 kg/m3, E 210 GPa, nu 0.3, yield 350 MPa | Placeholder |

Note on "clay with zero cohesion": undrained clay is modelled by cohesion alone (Tresca),
so with zero cohesion it would have no strength at all. A cohesionless bulk solid gets
its strength from friction, which is Mohr-Coulomb (or Drucker-Prager) with c = 0. T0
therefore uses the clay's density and stiffness with a friction angle.

Note on cost: with E = 10 MPa the P-wave speed is 87 m/s. A 5 mm cell at Courant number
0.3 gives a time step of about 17 microseconds, so 3 s of discharge is about 175,000
steps on about 13,000 cells. Lowering E (it barely affects granular flow) or mass
scaling are options if run time matters (WP 3.3p measures).

## Still open

- [ ] Do you have Abaqus or LS-DYNA access for reference runs? (Not needed for T0.)
- [ ] Licence: is GPL-3.0 acceptable (decision 0003)? (Needed before sharing the code.)
- [ ] Real material and tank data for T0, if you have a specific case.
