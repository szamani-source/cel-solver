# 0133 Initially fluidised dam break; the surface carries gas (WP F.8)

Date: 2026-11-24. Status: accepted; F.8, F.8b and F.8c not done (the case still misses the paper; split F.8d-F.8f) (follows decisions 0125, 0129, 0131, 0132).

## Context

F.8's done-when test: an initially fluidised dam break within 20 % of Roche et al. in
run-out and front speed, dry collapse with gas off unchanged, 2 and 4 processors equal
to serial to 1e-10 with gas mass to 1e-12, the cost with gas on against off (Sia's
rules, decision 0125 items 7-9; every test uses `ergunWenYu`).

The paper is Roche, Montserrat, Nino and Tamburrino (2008), J. Geophys. Res. 113,
B12203, doi 10.1029/2008JB005664 (HAL hal-00381402; decision 0125 and the F.0 table
called it 2010). Only its abstract could be read (IRD record fdi:010046767): HAL
(an automated-access block), AGU/Wiley, hal.science and the other copies refused
automated reading on 2026-11-24. The abstract gives: glass beads 60-90 micrometres,
columns initially fluidised and expanded 2.5-4.5 %, height/length 0.5-3; three phases
on the time scale sqrt(h0/g): a constant front velocity U = sqrt(2 g h0) from
t/sqrt(h0/g) = 1.5 to about 4, then a slower phase; the flow stops at about 6.5 with
run-out x/h0 = 5.5 to 6. Channel size, bead density, floor roughness and the gas supply
after release are not in it.

## Decision

1. **Case** (`tests/gas/collapse_cases.py`): a 0.2 x 0.2 m column (a = 1) against a
   smooth left wall, 10 mm cells (160 x 25), glass beads d 75 um, rho_p 2500, gas
   fraction 0.42 (0.40 loose packing expanded 3.5 %), plain Mohr-Coulomb (phi 30, c 0,
   psi 0, E 1 MPa), floor friction 30 deg, sealed floor, open top. It starts fully
   fluidised: zero effective stress and p = rho_b g (h0 - y), the gate gone at t = 0.
   Run-out read from the gate (x - x0)/h0. Assumed where the abstract is silent.
2. **The moving surface carries gas** (gasProperties `surfaceCarry yes`; default no,
   so earlier tests are untouched). Since F.4 a face between a gas cell and a void cell
   counted only grains crossing it, so when the surface sank inside its cell the
   vacated part had to fill with gas drawn through the grains (the Darcy drop over half
   a cell): the first run's surface cells sat at -1800 Pa (the Darcy drop for 0.2 m/s
   over 5 mm of 75 um beads is about 2200 Pa), which drained the column from the top and
   gave the grains friction. With `surfaceCarry` the surface face carries the
   Lagrangian volume flux like any other face, so the gas fills the vacated part (or
   leaves the newly filled part) from the void cell, conservatively (the void cell's
   balance takes the opposite). Faces between two void cells keep the old rule (gas flow
   there moves by its own velocity). Run-out 0.5 -> 1.3 h0. F.9 should decide whether
   it becomes the default after re-running the F.4-F.7 tests.
3. **The expanded bed is soft**: a fluidised bed has no grain contacts until it
   compacts back by its expansion; plain Mohr-Coulomb loads the skeleton on any
   compression (its constrained modulus, 1.35 MPa, is 5.6 times the gas's p0/e), so the
   grains took over the weight as soon as the column re-compressed. The case uses the
   compaction cap (decision 0113) normally consolidated from capP0 1 Pa with
   lambda* 0.0044 (3.5 % over ln 3000): a sealed bed then settles from the floor up at
   about U_mf, as a defluidising bed does. Run-out 1.3 -> 1.95 h0.

## Results

`tests/manual/test_fluidised_collapse.sh` (8 min on 4 cores, 2026-11-24; checkpointed,
results at /mnt/project-files/f8/):

| Check | Result | Limit | |
| --- | --- | --- | --- |
| Front speed, t/sqrt(h0/g) 1.5-4 | 0.589 m/s, -70 % of 1.981 | 20 % | FAIL |
| Run-out from the gate | 1.95 h0, -66 % of 5.75 | 20 % | FAIL |
| Stop | 6.3 sqrt(h0/g) | Roche 6.5 | |
| Gas mass (serial) | 1.2e-14 | 1e-12 | PASS |
| Dry collapse, gas off, against main (v1.1.0) | identical (0 difference) | 1e-12 | PASS |
| 2 processors (cut x = 0.1 m) | 1e-3 at 0.002 s, 5e-2 at the first write (0.02 s), blew up at 0.33 s | 1e-10 | FAIL |
| 4 processors (also y = 0.1 m) | 5e-2 at 0.02 s, blew up at 0.33 s; gas mass 3e-15 | 1e-10 | FAIL |
| Cost per step, gas on / off | 2.3x (223 s against 96 s, 9,945 steps, 4,000 cells) | F.9 target 1.5x | report |

Diagnostic runs (serial, same case): floor friction 0: 7.0 h0, front 1.28 m/s (the
flow reached the end wall at 1.6 m); internal friction 1 deg with the floor at 30 deg:
2.55 h0; the air kept on at U_mf through the reservoir floor: 2.6 h0, 0.58 m/s; without
the cap and without `surfaceCarry`: 0.5 h0, 0.05 m/s; dry: 0.65 h0, 0.116 m/s.

Lessons:

- What holds the flow back is the floor: the base defluidises first (sealed floor,
  pore gas leaving upward) and its friction on the floor then carries the load. Without
  floor friction the run-out overshoots. Roche's floor and gas supply after release
  are not in the abstract; the full paper is needed before tuning either.
- The parallel runs part from round-off within about 15 steps (1e-14 after 5 steps,
  1e-3 after 15), with or without `surfaceCarry`, and also with plain F.7 settings:
  every grain starts at the Mohr-Coulomb apex (zero effective stress), the worst case of
  the switch amplification decision 0132 saw after lift. A 2-processor run and a 5 mm
  serial run later blew up in the pore gas drag (`inertialGasDrag::slip`, square root
  of a negative number: the gauge pressure fell below -p0 after a void gas speed of
  150 m/s). Void gas speeds of 7-25 m/s appear in every run of this case (also before
  this package's changes): the void flow next to a spraying front is not robust.

## F.8b (2026-11-25)

Sia asked (through the coordinator) to carry on at once: find the paper's floor and gas
details, make the gas robust at the thin front, make 2/4 processors match.

1. **Paper**: still only abstracts. Roche (2008) and every citing paper tried
   (Montserrat et al. 2016 Granular Matter, Gueugneau et al. 2017 GRL and 2023 Bull.
   Volcanol., Breard et al. 2019 JGR) refuse automated reading or are paywalled. From
   their abstracts: 75 um beads; run-out x/x0 = alpha r^n, alpha growing with the degree of
   fluidisation and with bed expansion (up to more than twice non-fluidised); Breard
   et al. match the experiments with a frictional (Srivastava-Sundaresan) rheology. The
   floor and the gas supply after release stay assumptions; nothing was tuned.
2. **Crash fixed**: a nearly empty cell after the remap can hold more mass than its
   alpha (rho/alpha > rho_p); its filled-part gas fraction went negative, its
   permeability 0, the gas in it was shut in and its pressure diverged (to -p0, then the
   drag's square root failed). The filled-part gas fraction is now at least
   `eMinPacked` = 0.2, and the face gas densities in `faceCoefficients` use at least
   1e-3 p0. No run blows up now (serial, 2 and 4 processors to 1.3 s).
3. **Parallel fixed (the cap)**: serial and 2 processors parted by 1e-3 within 15 steps
   with the cap and not without it (gas and plain Mohr-Coulomb agreed to 1e-13). The
   cap's p_c of 0 means "not set"; the remap mixes it by mass (with a limiter) into
   loaded cells, so a cell could get p_c of 2e-5 Pa in one run and 1 Pa in the other,
   or 0 against 2765 Pa, and the next load then scaled its stress to almost nothing.
   `Cap::apply` now starts from p_c = max(stored p_c, P of the old stress, capP0), which
   is continuous in the stored value and, for a cell never loaded, the same as before.
   `tests/test_cap.sh`, `test_cap_stl.sh`, `test_settle_column.sh` and
   `test_constitutive_cpp.sh` still pass. It changes v1 results where the stress had got
   outside a diluted cap: the dry dam break with the cap moves its front by one cell.
   The v1 dry check therefore uses plain Mohr-Coulomb (version 1's collapse material).

Results (`tests/manual/test_fluidised_collapse.sh`, 2026-11-25, results in
/mnt/project-files/f8/run3, log test5.log):

| Check | F.8 | F.8b | Limit |
| --- | --- | --- | --- |
| Front speed | 0.589 m/s (-70 %) | 0.664 m/s (-67 %) | 20 % FAIL |
| Run-out | 1.95 h0 (-66 %) | 2.20 h0 (-62 %), stops at 7.8 sqrt(h0/g) | 20 % FAIL |
| Gas mass | 1.2e-14 | 1.0e-14 | 1e-12 PASS |
| v1 dry collapse against main | identical | identical (plain Mohr-Coulomb) | PASS |
| 2 processors | 5e-2 at 0.02 s, crashed 0.33 s | within 1e-10 to 0.36 s, 1.6e-8 at 1.3 s | 1e-10 FAIL |
| 4 processors | 5e-2 at 0.02 s, crashed 0.33 s | within 1e-10 to 0.46 s, 1.1e-8 at 1.3 s | 1e-10 FAIL |
| Cost gas on / off | 2.3x | 2.33x | 1.5x (F.9) |

The parallel difference now grows smoothly (1e-12 at 0.3 s, 1e-10 at 0.4-1 s, 1e-8 at
1.3 s, largest in Ug and sigma): round-off amplified by a spreading plastic flow, no
longer a jump. A restart from serial at 0.36 s stays within 1e-12 to 0.4 s.

What is left: the run-out (the defluidised base's floor friction; the paper's floor and
gas supply are unread), round-off growth over a long plastic flow against the 1e-10
limit, and the cost (F.9). Next: F.8c (`docs/work-packages.md`).

## F.8c (2026-11-25)

1. **Paper read** (Sia downloaded it; copies of every paper below in /mnt/project-files/refs).
   Roche et al. (2008): smooth Perspex apparatus, porous plate under the reservoir only ("no
   air flux is provided in the channel"), channel 3 m x 10 cm, h0 10-40 cm, x0 10-20 cm, d50
   75 um, rho_b 1450; Umf 7.5 mm/s (fully supported, not expanded), Umb 1.6-1.7 Umf
   (expanded 2.5-4.5 %). Run-out 5.5-6 h0 at Umb, 4.5-5 at Umf, 1.8-2.2 times dry.
2. **Case** set to it (`collapse_cases.py`): h0 0.4 m, x0 0.2 m (a = 2), 300 x 45 cells
   (3 m channel), floor 15 deg (Aravena et al. 2021), full support at t = 0, air kept on
   under the reservoir (patch floorRes, Darcy rate for full support).
3. **Remap overfill** (crash at 0.26 s): the remap clips alpha at 1 but keeps the mass, so a
   converging flow densified the packing with no stress (+17 % mean in F.8b, +20 % here, up to
   rho/alpha 2500) and squeezed the gas until the solve diverged. New opt-in
   materialProperties `overfillStress yes`: the overfill becomes an elastic volume strain
   ln(alpha) (mean stress -K ln alpha). Density now holds (+0.4 %). Also gasProperties
   `eMin` (the F.8b floor, default 0.2 unchanged). Neither changes v1 or earlier tests.
4. **Result**: run-out 1.25 h0 (90 % support: 1.20), front 0.39 m/s; gas mass 1e-14; v1 dry
   (plain MC) identical to main, 0.88 h0; 2 processors within 1e-10 to 0.06 s; 4 processors
   stopped on the project file limit.
5. **Diagnosis** (/mnt/project-files/f8/diag): support near the gate falls from 1 to 0 within
   4 ms and rings between -0.7 and +0.9 of the weight (gas stiffness p0/e ~ 240 kPa: 1 %
   dilation drops half the base weight); remap holes inside the bed (~0.8 % of its volume at
   20 ms); probes with 100x stiffer gas and permeability / 9 were worse. The v1 dry collapse
   itself is short: 0.85 h0 (MC 30 deg) against about 2.6 h0 in Roche's pressure-balanced
   runs (1.05 m, Breard et al. 2019) and about 2 h0 in Aravena's dry simulation; with mu(I)
   (tan 21, tan 33, I0 0.3) 1.2 h0, and with flowImpedance + fullBlend 3.5 h0.
6. **Published models** of these experiments: Gueugneau et al. (2017, depth-averaged Coulomb,
   advected basal pore pressure relaxing with a fitted D 0.004-0.02 m2/s, angles 24/19 deg),
   Aravena et al. (2021, 2D incompressible mu(I), advected-diffused p_f with fitted kappa
   0.015 against 0.13 theoretical, "reason unknown"), Breard, Dufek & Roche (2019, MFIX
   two-fluid, compressible gas and real permeability, kinetic theory plus the Princeton
   packing-dependent frictional pressure; delta 28, wall 12 deg, 5 mm, initial packing 0.58;
   the same pressure drop at release, recovered by compaction). Chupin & Dubois (2023,
   arXiv 2304.05083) is theory only.

Sia chose (2026-11-25): calibrate the dry collapse first, then route B (pore pressure as an
advected-diffused field with a fitted diffusivity, as Gueugneau and Aravena), then route A
(packing-dependent friction for the compressible model, as Breard). Split: F.8d (dry), F.8e
(B), F.8f (A).

## F.8d and F.8e (2026-11-25)

**F.8d, dry calibration.** The dry a = 2 column (gas off, no cap) with Gueugneau's angles
(phi 24 deg internal, floor 19 deg) and the mu(I) rule (mu2 tan 33, I0 0.3, d 75 um) plus
`flowImpedance` and `fullBlend` reaches 2.375 h0 from the gate (front at
1.15 m), -9 % of Roche's pressure-balanced 1.05 m (2.6 h0, Breard et al. 2019). Plain
Mohr-Coulomb (24/19 or 28/12) stops at 0.9-1.1 h0: its elastic face impedance over-damps a
fast collapse. `collapse_dry_reference.txt` was regenerated on main's binaries for this
material, so check B still proves version 2 leaves the dry solver unchanged.

**F.8e, route B.** New gasProperties `diffusivity` kappa (m2/s, default 0 = off, needs
`model poreGas`). With kappa > 0 the storage term becomes alpha k/(mu kappa) (so the pressure
diffuses at kappa, not at the rigid-bed k p0/(mu e)) and the squeezed-gas term is replaced by
upwind advection of p with the material that crosses each face in the remap: pore pressure
moves with the grains, with no compression or suction from the bulk volume change. This is
Iverson & Denlinger (2001) as used by Gueugneau (2017) and Aravena (2021). Gas mass is then
not conserved (pressure, not mass, is carried); the gas mass checks apply only to the
compressible route (kappa = 0). Test `tests/test_pore_diffusion.sh`: a step at a face of a
rigid bed matches erfc within 2 % at two kappa values, 2 and 4 processors within 1e-10.

Dam break (10 mm, serial), fitted kappa 0.015 (Aravena's): run-out 4.45 h0 (-23 % of 5.75),
front 2.09 m/s (-25 % of sqrt(2 g h0)), stops at 6.8 sqrt(h0/g) (Roche 6.5). With the
theoretical kappa k p0/(mu e) = 0.11: 3.30 h0, 1.18 m/s. For comparison the compressible route
with the same grains gives 3.15 h0 (with cap) and 2.57 h0 (without). Route B is the default
gas route of `collapse_cases.py` (KAPPA_FIT); kappa = 0 still runs the compressible case.

5 mm mesh (4 processors, kappa 0.015): run-out 5.12 h0 (-11 %, inside 20 %), but grains spray
ahead of the flow and made the mass-threshold front jump, so check A now takes the front as
the end of the continuous deposit (`deposit_history`). On that measure the front speed over
Roche's phase is 1.69 m/s at 5 mm (-40 %) and 1.48 m/s at 10 mm (-47 %; the 2.09 above
counted spray); 0.11: 1.18 m/s. Our flow is slower and lasts longer (stops at about 8
sqrt(h0/g) against 6.5). Run-out still depends on the mesh (4.45 to 5.12 h0). Plot
/mnt/project-files/f8/diag/routeB_front.png. Next F.8f.

## F.8f, first step (2026-11-25)

New materialProperties `princetonFr` (with `princetonR`, `princetonS`, `princetonMin`,
`princetonMax`, needs `rhoS`): the Princeton packing-dependent frictional pressure that
Breard et al. (2019) use, Pc = Fr (es - esMin)^r / (esMax - es)^s, es = rho/(alpha rhoS). It is
applied as a limit after the Lagrange step: an effective mean stress above Pc(es) is scaled
towards the apex (the grains cannot carry it at this packing, so the cell compacts); below
esMin the stress is zero. We did not take MFIX's form (Pc as the whole solid pressure, with
the dilation factor N) because our grains already have an elastic-plastic stress; the limit
keeps that and adds the packing dependence. Test `tests/test_packing_pressure.sh` (7 s): the
bound holds, a column at 0.58 compacts to 0.597 (static need 0.593: the drop from a loose
start overshoots), option off unchanged, 2 processors 1e-13.

Dam break, 10 mm, F.8d grains with Breard's law (Fr 0.05, r 2, s 5, 0.49-0.64) instead of
the cap: dry 2.47 h0 (Roche about 2.6); compressible gas (kappa 0) 3.50 h0, deposit front
1.54 m/s (was 3.15 h0 with the cap, 2.57 without); route B (kappa 0.015) 4.47 h0, 1.70 m/s
(4.45 without). The compressible run's base pore pressure under the reservoir falls from
5.6 to 0.9 kPa within 0.02 s and recovers to 3-4.6 kPa by 0.06-0.16 s, the behaviour Breard
reports (5.5 to 1 kPa, then about 4 kPa), but the pressure in the flowing part stays low
(0.2-1.6 kPa). Gap to Roche still -39 % for the compressible model.

Breard's own angles (28 deg internal, 12 deg floor) with our grains: dry 5.77 h0 (Roche about
2.6), compressible gas 7.00 h0. Our mu(I) + flowImpedance grains need the 19 deg floor to give
Roche's dry run-out (Breard's kinetic-theory solid dissipates differently), so the angles are
not taken over. Ratio fluidised/dry (Roche 1.8-2.2): compressible 1.42 (3.50/2.47), route B
1.81 (4.47/2.47).

## F.8g: route B as published (2026-11-25)

Diagnosis (Sia: use the papers to find why we are slower): with friction removed, the dry
column spreads at 1.43 sqrt(g h0), Roche's sqrt 2, so the numerics are not the limit. With
route B and the pressure kept (kappa 1e-6) the column interior stayed still and only the cell
at the free face was blown off (16 m/s): the carried p acted as a force, horizontally uniform
inside, all of its jump on the edge cell (the 5 mm spray). Gueugneau (2017) and Aravena (2021)
use the pore pressure only to lower the friction, tan phi = (1 - p/rho g h) tan delta.

Sia chose to implement that. New gasProperties `frictionOnly` (with `diffusivity`): sigma
stays the total stress, no gas force; each step lambda = min(p/P, 1) (P the total mean
pressure at the step start); the constitutive update sees sigma + lambda P_trial I (P_trial
the elastic trial pressure, so a supported cell keeps an elastic isotropic total stress) and
the wall friction is scaled by 1 - lambda. The route B advection now mixes the incoming
pressure by volume (bounded; the first runs diverged in fast thin flows); test_pore_diffusion
still passes. The case starts with the total stress -rho_b g (H0 - y) I.

Dam break 10 mm: fitted kappa 0.015 run-out 5.67 h0 (-1 %), front 1.04 sqrt(g h0) (-26 %;
Aravena's model 1.0-1.15); theoretical 0.11: 3.88 h0, 0.71; pressure kept: 1.25 (frictionless
1.43). `tests/manual/test_friction_only.sh` (coarse 20 mm): the supported flow is 25-29 %
behind the frictionless one: the carried p does not follow the dynamic rise of the total
pressure, so friction returns where P > p (the papers carry the ratio p/rho g h instead). 2
processors part from round-off at the Mohr-Coulomb apex, as the dry frictionless case does
(F.9).

5 mm (4 processors, kappa 0.015): front 1.18 sqrt(g h0) (-16 %, inside 20 %), at Roche's stop
time 6.5 sqrt(h0/g) the front is at 6.14 h0 (+7 %), but it keeps creeping and reaches the
channel end (7.0 h0) at 8.5. The fitted kappa depends on the mesh: the upwind advection of p
adds numerical diffusion that halves with the cell size (10 mm 5.67 h0, 5 mm > 7 h0).

## F.8h: support ratio carried (2026-11-25)

Sia approved carrying the support ratio and a sharper scheme. With frictionOnly the grains
now carry lambda = p/P (P the total mean pressure), as Gueugneau carries p/(rho g h): a cell
keeps lambda when its load changes, and material crossing a face brings its lambda,
reconstructed linearly to the face from the donor's gradient and bounded by the two cells'
values (second order, no new extrema). p = lambda P then diffuses with kappa. Two places
still used the full friction and are fixed: the mu(I) flow impedance (now on sigma +
lambda P I) and the wall Riemann solution of the stress predictor (mu (1 - lambda)).
`tests/test_friction_only.sh` (in the suite now): the fully supported coarse column moves
within 5.6 % of the frictionless one (25-29 % behind before).

Dam break 10 mm: theoretical kappa k p0/(mu e) = 0.11: run-out 5.02 h0 (-13 %), front 1.04
sqrt(g h0) (-26 %), stops at 5.6 sqrt(h0/g) (Roche 6.5); fitted 0.015 now overshoots (channel
end, front 1.26). So no fitted kappa is needed at 10 mm. 5 mm with the theoretical kappa:
front 1.19 sqrt(g h0) (-16 %, inside 20 %), at Roche's stop time 6.57 h0 (+14 %), then on to the
channel end (7.0 h0): run-out still depends on the mesh (10 mm 5.02 h0). Likely cause: the
void at p = 0 one half cell above the flow surface, so the drainage path of a thin flow scales
with the cell (Gueugneau relaxes p/(rho g h) over the depth). Refit on 5 mm (kappa 0.22) running.

## Consequences

- F.8 is not done. Its remaining work is split into F.8b (`docs/work-packages.md`):
  read Roche (2008) in full (setup, floor, gas supply, run table) and set the case to
  it; robust void flow and pore drag near a thin front (guard the drag against
  p < -p0, limit the void speed); parallel agreement from an apex start (the decision
  0132 issue, now the main blocker); then the 5 mm mesh.
- `surfaceCarry` and the soft-bed cap stay case options until F.8b shows they are right
  against the paper.
