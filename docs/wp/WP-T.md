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

- 2026-11-26: T.2b done. Temperature with immersed STL walls: capacity on the open
  volume, conduction on the open face area (STL walls adiabatic), T in the solid's
  small-cell mixing (by mass) and covered-cell ghosts. `tests/test_stl_heat.sh`
  (about 80 s, tutorials/stlHopper with a hot layer): energy 1.8e-14, T within
  20-150 C, uniform T 2e-13, heat out c T times discharged mass 4e-15. Left: T.7.
- 2026-11-26: T.8 done. `tests/test_thermal_restart_mpi.sh` (about 2.5 min): T.4
  column and T.6 bed restarted and on 2/4 processors within 3e-13; temperature
  1.0-1.1x a gas-on step (capacities and limits now computed once per step).
  Tutorial `tutorials/hotSiloDischarge` (T0 tank, 150 C layer under 20 C, cooled
  walls, Zehner-Schluender; README, plot.py; `tests/test_hot_silo_tutorial.sh`,
  22 s). Validation pack `tests/validation/thermal` (Allrun, about 5 min): T.4,
  T.5, T.6 on 1/2/4 processors. Manual: case setup, quick start, Validation
  (Temperature table), References. Lesson: oversubscribed MPI next to another
  run can take 40x longer (602 s against 14 s); time packs on a quiet machine.
  Left: T.7 (after the dam-break thread), T.2b (STL walls).
- 2026-11-26: T.6 done. Adapa et al. (2024) moving bed between plates:
  `tests/test_moving_bed_heat.sh` (about 50 s): local Nu within 0.45 % of the
  slug-flow solution with the gap, mean h 291.3 against 292.8 W/(m^2 K) (-0.52 %),
  inside their 225-350, energy 1e-14. `tools/thermal_reference.py` gains
  `slug_plates_nu_mean` and `slug_plates_h_mean` (eq. 9, checked against the length
  average of eq. 4 to 1e-6). The T limiter default is now MC (superbee steps smooth
  profiles: local Nu off by 20 %). Specific heat 1100 J/(kg K) assumed (not in the
  paper; h moves 1 % for 900-1300). Next: T.8 (T.7 waits for the dam-break thread).
- 2026-11-26: T.5 done, no solver change (the T.4 gas transport covers the void).
  `tests/test_void_heat.sh` (about 100 s): plates at a fixed wall temperature,
  Nu 7.457 against 7.5407 (-1.11 %, limit 3 %), steady to 1.6e-6, energy 3e-15.
  Next: T.6 (show Sia the mesh first).
- 2026-11-26: T.4 done. `model twoTemperature` (field `Tg`, Gunn exchange) and
  `equilibrium` with the gas; keywords `cGas`, `TgVoid`; poreGas accessors
  `massFlux()`, `qGas()`; celThermal links celFluid. `tests/test_gas_heat.sh`
  (about 30 s): Schumann gas 0.38 % and grains 0.26 % (limit 1 %), equilibrium
  front speed -0.02 %, uniform T under flow 1e-14, energy balance drift 2e-15,
  2 processors 2e-13. Next: T.5.
- 2026-11-26: T.3 done. thermalProperties `walls` (convective h, Tinf; heat flux q),
  `conductivity zehnerSchlunder` (`kParticle`, `kGas`, `rhoP`; field `kBed`).
  `tests/test_thermal_walls.sh` (about 25 s): slab against Carslaw-Jaeger 0.006 %
  (Bi 1) and 0.021 % (Bi 10), energy 2e-15, 2 processors identical, q A t to 6e-14,
  kBed against the reference tool 2e-16. Lesson: after changing a celThermal header
  rebuild celFoam too (Allwmake), or it crashes at exit. Next: T.4.
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
