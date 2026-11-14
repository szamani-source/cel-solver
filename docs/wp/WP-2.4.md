# WP 2.4: Remap in C++ (single material)

State: not started
Branch: `wp-2.4-remap`

## Goal

`eulerianSolid::remap(dt)` returns the Lagrange-step state to the fixed mesh:
conservative, second order, monotone, per-field limiters, on any `fvMesh`, serial and
parallel.

## Read first

- `docs/decisions/0014-lagrange-remap-contract.md` (exactly what remap receives and must
  return), `0013` item 5 (limiters), `0008` item 2
- `prototype/lr1d.py` and `prototype/lr2d.py`: the 1D swept-region MUSCL remap
- `src/celEulerian/eulerianSolid/`, `tests/constitutive/wrapperCheck/` (pattern for a
  small OpenFOAM test program)
- `CLAUDE.md`: run `tools/setup_openfoam.sh` first

## Build

In `src/celEulerian/eulerianSolid/eulerianSolidRemap.C` (new file; move the stub there):

1. Input: Lagrangian fields and `phiL_`. Lagrangian volume `V_L = V + dt * sum(phiL)`.
2. Per face the swept volume `dV = dt * phiL`; donor cell = the cell it leaves.
3. Limited linear reconstruction in the donor cell, evaluated at the centroid of the
   swept volume (face centre minus half the face displacement): density for the mass
   flux `dm`; velocity, stress and history per unit mass, fluxed with `dm`
   (mass-weighted), so `sum(m q)` is conserved exactly.
4. Limiter selectable per field in `system/fvSchemes` or `materialProperties`
   (thread's choice, record it), defaults MC for density, velocity and stress and
   superbee for history. On a general mesh use a cell-gradient limiter that reduces to
   these in 1D on a uniform mesh; say in the decision note what was used.
5. New state on the fixed mesh: `rho = (m - sum dm)/V`, `q = (m q - sum dm q_f)/m_new`.
   Cells with no flux through any face are left bit-identical.
6. Boundary faces: zero flux through walls, symmetry and empty; inflow and outflow patches
   out of scope unless trivial; processor faces must give identical fluxes on both sides.
7. `tests/eulerian/remapCheck/`: a test program that builds an `eulerianSolid`, sets the
   fields by hand through the `...Ref()` accessors, sets `phiL` from a prescribed velocity
   field and the density consistently (`rho_L = rho V / V_L`), and calls `remap` many
   times. This tests remap alone, without WP 2.3.

`tests/test_remap.sh`.

## Done when

`tests/test_remap.sh` passes in CI (remapCheck with prescribed motion):

- Uniform translation along x on a 3D mesh reproduces the 1D prototype remap of the
  same profile and Courant number to 1e-10 (the 3D code reduces to 1D).
- Smooth profile translated diagonally in 2D and along the body diagonal in 3D for one
  box crossing: order 1.7 or better over 3 meshes; state the largest stable material
  Courant number found and assert at the default.
- Square profile, diagonal translation: no new extrema above 1e-10 for density, velocity
  and stress; history edge width (10 to 90 percent) at most 4 cells after 50 cells of travel.
- Solid-body rotation of a slotted disc of history, one turn: total preserved to
  round-off, no new extrema above 1 percent, shape error reported.
- Mass, momentum and `sum(m q)` for every field conserved to 1e-12 relative in all cases.
- A uniform state stays uniform to round-off under a non-uniform divergence-free flux
  and under a compressing flux followed by its reverse.
- The diagonal case on 2 and 4 processors matches serial to 1e-12.
- Existing tests still pass unchanged.

Note: without corner transport an unsplit donor-cell remap is stable only for a summed
Courant number below about 1; the acoustic time step (about 0.3 c) keeps material
Courant numbers far below that except for supersonic motion. Record what was measured.

## Out of scope

- `lagrangeStep` (WP 2.3, running in parallel in `eulerianSolidLagrange.C`; both touch
  `Make/files` and `docs/STATUS.md`). Combined checks: WP 2.4b.
- Volume fractions and interface reconstruction (WP 3.1), the `apex` flag (not stored,
  decision 0010).

## Log

(none yet)
