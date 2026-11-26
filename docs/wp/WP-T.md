# WP T.0 to T.8: temperature (energy equation)

State: in progress
Branch: `claude/project-thread-8aqh4z` (thread "Add temperature to the solver")

## Goal

An energy equation for the grains and the gas: conduction, heat carried by the
moving solid and the gas, grain-gas exchange (Gunn), wall heating and cooling,
restart and MPI. Plan approved by Sia 2026-11-26 (decision 0135).

## Read first

- `docs/decisions/0135-temperature-plan.md`, `docs/spec/thermal-references.md`
- Gas side: decisions 0126 (fractions), 0127 (pore pressure), 0129 (void, outlets),
  0131 (void gas flow); `src/celFluid/poreGas/poreGas.H`
- Remap of a mass-weighted history field: the cap's p_c (decision 0113) in
  `src/celEulerian/eulerianSolid/eulerianSolidRemap.C`

## Build

- `src/celThermal` (library), `constant/thermalProperties`, fields `T` and `Tg`.
- One test per package (table "Temperature (T series)" in `docs/work-packages.md`).

## Done when

- Each package's test passes; T.8 passes restart and 2/4 processor checks and the
  manual check.

## Out of scope

- Thermal stress, temperature-dependent friction, radiation beyond Zehner-Schluender,
  drying and phase change, reactions (decision 0135).

## Log

- 2026-11-26: T.2 done. `eulerianSolid::addPassenger` and `storeRemapMass`;
  heatTransfer `beforeRemap`/`afterRemap`, `limiter` (default superbee), `voidT`;
  celFoam prints "Heat carried in". `tests/test_carried_temperature.sh` (about 30 s):
  block carried 50 cells, energy 2e-14, front in place (1.1011 against 1.1 m) and
  2 cells wide, 4 processors 2e-14; box with inlet (353 K) and outlet, balance
  8e-15, fed material at the inlet temperature exactly; one temperature exact;
  carried heat = c T (fed - out) to 1e-13; 2 processors 2e-15 until the solid parts
  (decision 0134). Split: T.2b immersed STL walls with temperature. Next: T.3.
- 2026-11-26: T.1 done. `src/celThermal` (`heatTransfer`), `constant/thermalProperties`
  (`model equilibrium`, `c`, `conductivity constant`, `k`), field `T`, celFoam calls it
  after the gas and prints thermal energy, heat in and their balance at each write.
  `tests/test_conduction.sh` (about 15 s): erfc 0.012 % (limit 1 %), heat in -0.024 %,
  constant flux 0.006 %, linear temperature on 27 deg cells exactly steady, energy
  7e-15 adiabatic and 3e-15 with a heated face, 2/4 processors 6e-15, off = byte for
  byte. Lessons in decision 0135 (processor patch values, patch delta). Manual: new
  thermalProperties section in case-setup.html. Next: T.2 (temperature carried by
  the remap).
- 2026-11-26: T.0 done. `tools/thermal_reference.py` (pure Python) and
  `tests/test_thermal_reference.sh` (about 20 s): erfc and the cooled slab against
  finite differences and Incropera Table 5.1, Zehner-Schluender limits and bounds, Gunn
  limits, Schumann against an independent finite-difference solution (2e-6) and the
  equilibrium front, Graetz plates 7.5407 against 7.541, slug flow with a wall gap
  against marching (1e-4). Original Gunn, Zehner-Schluender and Kuipers papers are not
  readable from the container; formulas taken from three open papers that agree.
  Sullivan and Sabersky unreadable: T.6 now uses Adapa et al. (2024). Next: T.1.
