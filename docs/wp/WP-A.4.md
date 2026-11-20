# WP A.4: Sand footing with psi 30

State: done (2026-11-19), decision 0095. Needs A.3b. Thread "Fix sand footing bearing
capacity". Split: the remaining N_gamma excess goes to A.5.

## Goal

Find why the sand footing (decision 0057) reads high and why associated flow blows up.

## Done when

`tests/test_sand_dilation.sh`: phi 30, psi 30 footing at 10 mm pushed to 0.06 B runs
through, q / q_ref below 3 and steady over the second half, mass to 1e-10, density at
most 1.01 rho0.

## Log

- 2026-11-19: domain width ruled out. The blow-up is a partly filled surface cell that
  cannot dilate: `voidPlasticFlow` fixes it. The excess is a numerical rate effect
  (rho c v): psi 30 reads 20 / 16.6 / 13.9 / 11.9 as the push halves three times, psi 0
  10.9 / 9.3. Manual N_gamma check still fails; next A.5 (rate-free Mohr-Coulomb flow
  impedance, psi 30 locking past 0.06 B).
