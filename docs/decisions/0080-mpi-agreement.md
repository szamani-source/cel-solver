# 0080 What "1 and N cores agree" means (WP 6.1)

Date: 2026-11-18. Status: accepted.

## Context

WP 6.1's done-when is "1 and N cores agree" for runs with Lagrangian bodies. The Eulerian
packages already checked parallel runs cell by cell (2.3, 2.4, 2.4b, 3.2a: 1e-13 to
1e-15), but the collapse could not (3.2b, 3.4b: round-off flips the shear cap, the void
rules or the bisection and the runs part). Bodies add more such switches: cut-cell
mixing, hand-overs, contact on and off.

First runs (2026-11-18): the footing push on 2 processors matches serial to 1e-14 for
221 steps, then parts by 1e-11 to 1e-7 in two steps where the footing enters the top cell
row. A serial run with the footing speed changed by 4e-14 (relative) parts by the same
amount at the same step (8e-7 against 3e-7). So the jump is the scheme amplifying
round-off at a switch, not a parallel error. Sums over processors (body force, node
forces of deformable bodies, mass) are added in a different order, so parallel runs
always carry round-off differences.

## Decision

A parallel check has two parts (`tests/test_mpi_bodies.sh`):

1. **Exact phase.** The first 100 steps of the history (mass, momentum, body pose,
   velocity and force, wall forces, largest stresses) match serial to 1e-12 of each
   column's largest value (or 10 times the difference of a perturbed serial run where
   that is larger: a stress-free fill sits at yield and its cap counts follow round-off
   from step 2). A missing reduction or a wrong value on a processor boundary
   shows here at once.
2. **Outcome.** The case result is within a tolerance of serial, and a serial run with
   maxCo changed by 1 part in 1e12 is printed beside it as the round-off band:
   - free block dropped onto T0 (2 processors through the block): force on the block at
     rest within 2 %, height within 0.2 mm;
   - steel hex8 footing pushed into T0 (2 processors through it): mean force within 5 %
     (the WP 5.3a tolerance; the perturbed serial run alone moves it 2.4 %);
   - T0 tank with STL walls, floor and outlet (4 processors, 2 x 2, splits through the
     outlet and the floor): discharge rate within 5 % (the WP 4.5 tolerance).
   Mass to 1e-10 and no hand-over loss on N processors; the tank also no leak while
   closed and tank plus discharged constant.

## Consequences

- Parallel results with bodies are the same physics as serial, not the same digits.
  Compare parallel and serial runs as bands (as for the silo rate, decision 0046).
- No code changed for MPI: bodies are whole on every processor (STL, hex8 elements),
  forces are summed with `returnReduce` / `listCombineReduce`, cut cells are local.
- Small cells mix only with neighbours on their own processor (decision 0051 mixes
  through faces of the cell; processor faces are not used). No case showed an unmixed
  small cell at a processor boundary (`nUnmixed` 0); a body lying along a processor
  boundary could leave one unmixed, which the solver reports.
