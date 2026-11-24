# 0133 Initially fluidised dam break; the surface carries gas (WP F.8)

Date: 2026-11-24. Status: accepted, package not done (follows decisions 0125, 0129, 0131, 0132).

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

## Consequences

- F.8 is not done. Its remaining work is split into F.8b (`docs/work-packages.md`):
  read Roche (2008) in full (setup, floor, gas supply, run table) and set the case to
  it; robust void flow and pore drag near a thin front (guard the drag against
  p < -p0, limit the void speed); parallel agreement from an apex start (the decision
  0132 issue, now the main blocker); then the 5 mm mesh.
- `surfaceCarry` and the soft-bed cap stay case options until F.8b shows they are right
  against the paper.
