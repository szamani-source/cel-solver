# Target problems, units and material data (WP 0.3)

Status: DRAFT. Defaults below were chosen by Claude; Sia to confirm or change.

## Units

SI throughout: m, kg, s, Pa. Stress positive in tension. Gravity -9.81 m/s^2 in z.

## Target problems (proposed defaults)

| # | Problem | Why | Lagrangian body | Eulerian material | First version |
| --- | --- | --- | --- | --- | --- |
| T1 | Rigid strip footing on undrained clay | Analytic answer (Prandtl, Nc = 5.14) | Rigid | Tresca clay | Plane strain |
| T2 | Cone or pile penetration in clay | The classic CEL use; Abaqus and MPM results published | Rigid | Tresca, then Mohr-Coulomb | 3D quarter model |
| T3 | Projectile into soil | Exercises deformable bodies and high rates | Deformable (J2 steel) | Mohr-Coulomb or Drucker-Prager | 3D |

## Material data (placeholders)

| Material | Model | Parameters |
| --- | --- | --- |
| Undrained clay | Tresca (J2 with su) | rho 1800 kg/m3, E 10 MPa (E/su about 500), nu 0.49, su 20 kPa |
| Sand | Mohr-Coulomb | rho 1900 kg/m3, E 50 MPa, nu 0.3, phi 35 deg, c 1 kPa, psi 5 deg |
| Steel | J2, later Johnson-Cook | rho 7850 kg/m3, E 210 GPa, nu 0.3, yield 350 MPa |

## Open questions for Sia

- [ ] Are T1 to T3 the right problems, or do you have specific cases (geometry, data)?
- [ ] 2D (plane strain, axisymmetric) first, or straight to 3D?
- [ ] Do you have Abaqus or LS-DYNA access for reference runs?
- [ ] Licence: is GPL-3.0 acceptable (decision 0003)?
