# Evolution Equation Map

Use `docs/raw/derivation.md` for equation authority. This page maps equation concepts to code and tests without repeating the derivation.

## Conservative Equations

Concept:
- Evolution is written for conservative variables rather than primitive variables.
- The derivation summarizes density, energy, momentum, and magnetic-field conservative equations.

Code map:
- Conservative struct: `GRHayL/include/ghl.h`
- Primitive-to-conservative conversion: `GRHayL/Con2Prim/compute_conservs.c`
- Conversion plus stress-energy: `GRHayL/Con2Prim/compute_conservs_and_Tmunu.c`
- Conservative limiting: `GRHayL/Con2Prim/apply_conservative_limits.c`
- Undensitization for recovery: `GRHayL/Con2Prim/undensitize_conservatives.c`
- Con2Prim helper route: [limits and conversions](../gems/con2prim/limits-and-conversions.md)

Tests:
- `Unit_Tests/unit_test_compute_conservs_and_Tmunu.c`
- `Unit_Tests/unit_test_apply_conservative_limits.c`
- `Unit_Tests/unit_test_ET_Legacy_conservs.c`
- `Unit_Tests/data_gen/unit_test_data_ET_Legacy_conservs.c`

## Stress-Energy Tensor

Concept:
- `docs/raw/derivation.md` derives the electromagnetic contribution and combined stress-energy tensor with GRHayL magnetic rescaling.
- Stress-energy feeds conservative energy, momentum, fluxes, and source terms.

Code map:
- Public struct: `ghl_stress_energy` in `GRHayL/include/ghl.h`
- Small-b and magnetic scalar: `GRHayL/GRHayL_Core/compute_smallb_and_b2.c`
- Tensor helpers: `GRHayL/GRHayL_Core/compute_TDNmunu.c`, `GRHayL/GRHayL_Core/compute_TUPmunu.c`
- Conservative and tensor combined path: `GRHayL/Con2Prim/compute_conservs_and_Tmunu.c`

Tests:
- `Unit_Tests/pert_test_fail_stress_energy.c`
- `Unit_Tests/unit_test_compute_conservs_and_Tmunu.c`
- `Unit_Tests/data_gen/unit_test_data_grhayl_core_test_suite.c`

## Flux Terms

Concept:
- Flux terms move conservative quantities through cell faces.
- HLLE flux paths are split by direction, EOS family, and entropy support.
- Reconstruction provides face values before flux evaluation.

Code map:
- Public flux API: `GRHayL/include/ghl_flux_source.h`
- Characteristic speeds: `GRHayL/Flux_Source/ghl_calculate_characteristic_speed_dirn0.c`, `GRHayL/Flux_Source/ghl_calculate_characteristic_speed_dirn1.c`, `GRHayL/Flux_Source/ghl_calculate_characteristic_speed_dirn2.c`; route caller contract through [Flux_Source characteristic speeds](../gems/flux-source/characteristic-speeds-contract.md)
- Hybrid, hybrid entropy, tabulated, and tabulated entropy HLLE variants: route
  direction/EOS variant mapping through [Flux_Source HLLE flux variant matrix](../gems/flux-source/hlle-flux-variant-matrix.md)
- Reconstruction sources: `GRHayL/Reconstruction/`

Availability boundary:

- Family-specific direct HLLE functions, legacy `void` names and their
  `_checked` counterparts, are the only callable public HLLE routes. Deprecated
  generic `ghl_calculate_HLLE_fluxes_dirn0/1/2` pointers remain declared in
  `GRHayL/include/ghl_eos_functions.h` and exported, but EOS initialization
  never assigns them.
- No-HDF5 builds retain the direct tabulated HLLE definitions; their
  table-dependent tests and generators remain excluded.
- Reconstruction is built public caller API, but no production source caller is
  visible in this repository. Tests demonstrate family kernels with documented
  assertion/generator limits; they do not prove an evolution integration.

Tests:
- `Unit_Tests/unit_test_hybrid_flux.c`
- `Unit_Tests/unit_test_tabulated_flux.c`
- `Unit_Tests/unit_test_ET_Legacy_flux_source.c`
- `Unit_Tests/unit_test_tabulated_eos_compose.c`
- `Unit_Tests/data_gen/unit_test_data_hybrid_flux.c`
- `Unit_Tests/data_gen/unit_test_data_tabulated_flux.c`

Induction HLL tests (`unit_test_HLL_flux.c`, `unit_test_ET_Legacy_HLL_flux.c`)
route through the Induction section below, not Flux_Source HLLE coverage.

Existing flux tests cover wave-bound and callback-error contracts. `make tests`
compilation or workflow selection alone is not runtime proof.

## Source Terms

Concept:
- Source terms include spacetime-derivative contributions to hydrodynamic evolution.
- Callers provide metric derivatives in the form expected by the source routine.

Code map:
- Public source API: `GRHayL/include/ghl_flux_source.h`
- Source implementation: `GRHayL/Flux_Source/ghl_calculate_source_terms.c`;
  route caller-owned metric derivative and extrinsic-curvature contract through
  [Flux_Source source-term contract](../gems/flux-source/source-terms-contract.md)
- Stress-energy inputs: `GRHayL/GRHayL_Core/compute_TDNmunu.c`, `GRHayL/GRHayL_Core/compute_TUPmunu.c`
- NRPy 2 source generator: `GRHayL/Flux_Source/generate_flux_source.py`;
  route generation through `GRHayL/Flux_Source/generate_flux_source.sh` and
  [Flux_Source generated NRPy boundary](../gems/flux-source/generated-nrpy-boundary.md)

Tests:
- `Unit_Tests/unit_test_ET_Legacy_flux_source.c`
- `Unit_Tests/data_gen/unit_test_data_ET_Legacy_flux_source.c`
- `Unit_Tests/pert_test_fail_conservatives.c`

## Radiation M1 Moments And Neutrino Transport

Concept:
- The Radiation gem evolves grey, one-group, three-species neutrino number and
  energy/momentum moments alongside the GRMHD conservative state; it is a
  library-level gem, not a downstream host evolution route.
- Equation authority stays in the Radiation leaves and Doxygen source: route
  the 3+1 radiation equations through
  [M1 3+1 radiation equations](../gems/radiation-m1/m1-3plus1-radiation-equations.md),
  the four-dimensional closure through
  [M1 four-dimensional Minerbo closure](../gems/radiation-m1/m1-four-dimensional-minerbo-closure.md),
  neutrino number-current and transport terms through
  [M1 number current and transport](../gems/radiation-m1/m1-number-current-and-transport.md),
  comoving moments and source projection through
  [M1 comoving moments and source projection](../gems/radiation-m1/m1-comoving-moments-and-source-projection.md),
  matter/lepton exchange through
  [neutrino lepton and matter exchange](../gems/radiation-m1/neutrino-lepton-and-matter-exchange.md),
  and neutrino source equations through
  [M1 neutrino source equations](../gems/radiation-m1/m1-neutrino-source-equations.md).

Code map:
- Public M1 variables, parameters, and operation declarations:
  `GRHayL/include/ghl_m1.h`; neutrino rate-provider declarations:
  `GRHayL/include/ghl_neutrino_rate_provider.h`; aggregate include:
  `GRHayL/include/ghl_radiation.h`. Installed-header ownership routing:
  [Public API Map](../public-api-map.md) and
  [Radiation M1 API/build boundary](../gems/radiation-m1/api-build-boundary.md).
- Initialization and helper-API routing starts at the
  [M1 public API quick reference](../gems/radiation-m1/m1-public-api-quick-reference.md)
  leaf, which indexes the installed `ghl_m1.h` initialization, solver/closure
  control, and shared helper surface; shared kernels live in
  `GRHayL/Radiation/` with neutrino operations in
  `GRHayL/Radiation/Neutrinos/`.
- Face transport: canonical four-point blended Rusanov
  (`GRHayL/Radiation/ghl_m1_four_point_blended_rusanov.c`) and the shared
  component-wise Rusanov helper declared in `ghl_m1.h`; flux details route
  through [M1 finite volume and face flux](../gems/radiation-m1/m1-finite-volume-and-face-flux.md)
  and [M1 four-point blended Rusanov](../gems/radiation-m1/m1-four-point-blended-rusanov.md).
- Local coupling: explicit/frozen-rate/pair source updates and the implicit
  solve route through
  [M1 source update branches and rollback](../gems/radiation-m1/m1-source-update-branches-and-rollback.md)
  and [M1 implicit Newton and failure policy](../gems/radiation-m1/m1-implicit-newton-and-failure-policy.md).
- Radiation stress-energy inputs to any host `Tmunu` assembly route through
  `GRHayL/Radiation/ghl_m1_stress_energy.c`; the host-stage boundary routes
  through [M1 host stage and volume-weighted integration](../gems/radiation-m1/m1-host-stage-and-volume-weighted-integration.md).

Tests:
- `Unit_Tests/unit_test_m1_closure_fallback.c`
- `Unit_Tests/unit_test_m1_diffusion_flux.c`
- `Unit_Tests/unit_test_m1_error_handling.c`
- `Unit_Tests/unit_test_m1_fd_jacobian.c`
- `Unit_Tests/unit_test_m1_neutrino_rusanov_flux.c`
- `Unit_Tests/unit_test_m1_neutrino_seeded_invariants.c`
- `Unit_Tests/unit_test_m1_neutrino_source_update.c`
- `Unit_Tests/unit_test_m1_rate_provider.c`
- `Unit_Tests/unit_test_m1_thcm1_blended_rusanov.c`
- `Unit_Tests/unit_test_rusanov_flux.c`
- Fixture, runner, and coverage routes: [M1 tests and fixtures](../gems/radiation-m1/tests-and-fixtures.md),
  [Test Map](../test-map.md), and the [M1 unit-test guide](../../docs/raw/Radiation_unit_tests.md).

Availability boundary:

- `configure` builds the M1 targets with the rest of the unit tests; the
  ordinary `.github/run_tests.sh` invokes them with pinned TestData replay by
  default; `M1_FIXTURE_DIR` selects a supplied package instead of downloads.
  The gem is library-level: no host evolution,
  Cactus build, or physical validation is established by these tests.
- Radiation is not part of the hydrodynamic Flux_Source HLLE family; do not
  count its Rusanov helpers as Flux_Source coverage.

## Induction and Vector Potential

Concept:
- Magnetic evolution uses vector potential and densitized scalar potential paths documented in `docs/raw/Induction.dox`.
- The induction gem separates magnetic HLL flux terms from gauge/scalar-potential RHS terms.
- Keep equation authority in `docs/raw/Induction.dox` and `docs/raw/derivation.md`; use KB pages only for routing.

Code map:
- Public induction API: `GRHayL/include/ghl_induction.h`
- Magnetic HLL flux: `GRHayL/Induction/HLL_flux_with_B.c`, `GRHayL/Induction/HLL_flux_with_Btilde.c`; route caller contract through [Induction HLL flux contract](../gems/induction/hll-flux-contract.md)
- Interpolation: `GRHayL/Induction/Interpolators/`
- Scalar-potential RHS: `GRHayL/Induction/calculate_phitilde_rhs.c`; route gauge/scalar RHS contract through [Induction gauge RHS contract](../gems/induction/gauge-rhs-contract.md)
- Characteristic speed dependency:
  [Flux_Source characteristic speeds](../gems/flux-source/characteristic-speeds-contract.md)

Tests:
- `Unit_Tests/unit_test_induction_ccc_ADM.c`
- `Unit_Tests/unit_test_induction_ccc_BSSN.c`
- `Unit_Tests/unit_test_induction_vvv_ADM.c`
- `Unit_Tests/unit_test_ET_Legacy_induction_gauge_rhs.c`
- `Unit_Tests/unit_test_HLL_flux.c`
- `Unit_Tests/unit_test_ET_Legacy_HLL_flux.c`
- `Unit_Tests/data_gen/unit_test_data_HLL_flux.c`
- `Unit_Tests/compute_A_flux_with_B.c`
- `Unit_Tests/compute_A_flux_with_Btilde.c`
- `Unit_Tests/data_gen/unit_test_data_induction_interpolation.c`
- Fixture and run routes: [Induction tests and fixtures](../gems/induction/tests-and-fixtures.md),
  [Induction verification workflows](../gems/induction/verification-workflows.md)

## Conservative-to-Primitive Recovery

Concept:
- Conservative evolution requires recovering primitive variables after updates.
- Recovery is numerical and EOS-dependent.

Code map:
- Public recovery API: `GRHayL/include/ghl_con2prim.h`
- Multi-method dispatch: `GRHayL/Con2Prim/con2prim_multi_method.c`
- Hybrid methods: `GRHayL/Con2Prim/Hybrid/`
- Tabulated methods: `GRHayL/Con2Prim/Tabulated/`
- Guess and bounds: `GRHayL/Con2Prim/guess_primitives.c`, `GRHayL/Con2Prim/enforce_primitive_limits_and_compute_u0.c`
- Diagnostics: `GRHayL/Con2Prim/initialize_diagnostics.c`
- KB routes: [solver matrix](../gems/con2prim/solver-matrix.md),
  [recovery flow](../gems/con2prim/recovery-flow.md),
  [limits and conversions](../gems/con2prim/limits-and-conversions.md)

Tests:
- `Unit_Tests/unit_test_con2prim_multi_method_hybrid.c`
- `Unit_Tests/unit_test_con2prim_tabulated.c`
- `Unit_Tests/unit_test_con2prim_debug.c`
- `Unit_Tests/unit_test_hybrid_failure.c`
- `Unit_Tests/data_gen/unit_test_data_con2prim_multi_method_hybrid.c`
- Fixture route: [Con2Prim tests and fixtures](../gems/con2prim/tests-and-fixtures.md)

Recovery floor note: `apply_conservative_limits` sets and accumulates
`diagnostics->tau_fix` when raising tau to `tau_atm`; family-dependent
`tau_atm` setup is routed through [EOS dispatch](../core/eos-dispatch-contract.md).
