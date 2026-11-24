# 0129 Pore gas with a free surface, void and outlets (WP F.4)

Date: 2026-11-24. Status: accepted (follows decisions 0125 item 4, 8, 9 and 0128).

## Context

F.3 coupled the gas and the solid both ways but stopped with an error for a
`freeSurface` solid: gas in void cells had no mobility, so material moving in would
trap and compress it. Decision 0125 item 4 asks for void cells at the boundary
pressure, partly filled cells with k from the filled part, and the free surface and
outlet faces carrying the void pressure. Sia asked that every test use `ergunWenYu`,
that 2 and 4 processors match serial to 1e-10 with gas mass to 1e-12, with cuts through
the bed, the free surface and the outlet, and that the cost with gas on be reported
against the 1.5x target.

## Decision

1. **Void cells hold the void pressure.** With a `freeSurface` solid, a cell whose
   alpha is below `alphaVoid` (gasProperties, default 1e-4, the v1 void scale of
   `fullBlendVoid` and the nearly-full rule) is void for the gas: its p is 0 gauge (the
   ambient p0), fixed in the solve (no matrix coupling to it, so the matrix stays
   symmetric across processors too). The atmosphere gives or takes what its balance
   lacks; the bookkeeping counts that as inflow (`voidIn`, part of `netIn`), so gas mass
   less inflow stays constant to round-off. Gas flow inside the void is F.6.
2. **The free surface is a pressure face.** A face between a gas cell and a void cell
   carries the gas cell's mobility over the distance from its centre to the face, like a
   fixedValue boundary (the venting check: within 0.3 % of Terzaghi). The solid sees the
   void pressure on that face (`poreGas::pFace()`, passed with
   `eulerianSolid::setPorePressure(p, pf)`), so the uplift of a bed is p at its base less
   0 at its surface, not less half the top cell's p.
3. **The free surface is open.** At a face between a gas cell and a void cell no gas is
   squeezed by the face's motion: only the grains crossing it displace their volume of
   gas (G = -dm/rhoP). Counting the mixture volume (dt phiL - dm/rhoP, as on other faces)
   let a creeping surface cell pump void gas into the bed: the vented column read 14 %
   low (stiff bed) and 31 % low (soft bed). v1's partly filled surface cell creeps when
   the surface settles by more than about 1e-4 of a cell, with or without gas (known
   issue).
4. **Walls carry no gas.** Where material pulls away from a wall (v1 opens void there)
   no gas crosses the wall; the cell's gas fills the gap. Only faces open to the gas
   outside carry it with the material: open outlets and inlets and fixedValue patches of
   p (out at the cell's pressure, in at the boundary's). Before this rule round-off
   wall-normal velocities brought 1e-9 kg in through the walls of the outlet box before
   its outlet opened.
5. **Outlets.** A solid outlet (materialProperties `outlets`) must be fixedValue in
   0/p (without 0/p it is fixedValue 0); it is sealed for the gas until it opens. The
   remap keeps the mass through outlets for the gas (`remapMass` on non-coupled
   patches; `dmPatch` for `remapField` unchanged), so gas leaves with the material.
6. **Partly filled cells** keep k of their filled part (F.1) and store the gas of their
   whole open volume.

## Consequences

- `tests/test_gas_venting.sh` (about 2 min) passes: column of fine powder (d 0.1 mm, e
  0.40, E 100 MPa) 0.5 m, 21 cells, 7 void cells above, gas at 100 Pa over p0: venting
  within -0.28 / -0.14 / -0.04 / -0.05 % of Terzaghi at T = 0.05 / 0.1 / 0.2 / 0.5
  (limit 3 %), void cells at p0 exactly, gas mass 7e-16. The soft bed of F.3 (E 1 MPa)
  reads +2.8 % at T = 0.05 falling to +0.6 % at 0.5 (report: its surface cell creeps).
  Outlet box (Mohr-Coulomb, 1 mm grains): no gas crosses the boundaries before the
  outlet opens (exactly 0), the outlet holds p0, gas mass 2e-15, discharge with gas
  within 0.1 % of the dry run (coarse grains). 2 and 4 processors, cuts across the
  column in the bed and at the surface, and through the outlet and the bed of an
  elastic box and the Mohr-Coulomb box (material leaving): p within 3e-11, U within 1e-10, gas mass 1e-15
  on every count.
- Cost per step, gas on against off: 2.6 to 2.7x on the 28-cell column (2 solves per step,
  per-step overhead), 1.9x on the 480-cell outlet box, 2.2 to 2.3x at 2.5 mm (7,680
  cells, 4.6 solves per step: Ergun inertia at 1 mm grains and a stiff diffusion,
  c dt/dx^2 about 100). Above the 1.5x target of F.9. A profile (callgrind, 1,920
  cells) puts the gas at half the run: PCG 21 %, face coefficients 12 % (four drag-law
  calls per face), two matrix assemblies per iteration 9 %, non-orthogonal correction
  on an orthogonal mesh 4.5 %. F.9 should merge the secant and tangent mobility into one
  call, skip the correction on orthogonal meshes and assemble once per iteration; GAMG
  (19 iterations per solve against 100 for PCG) was slower here (setup each solve).
- v1 unchanged: without gas nothing in the remap or the momentum changes
  (`test_void`, `test_outlet` pass); F.1 to F.3 tests pass unchanged (no free surface).
- Not yet: gas flow in the void (F.6), inlets beyond carrying the gas in at the boundary
  pressure (F.5), restart of the void bookkeeping (`netIn` restarts at 0, F.9).
