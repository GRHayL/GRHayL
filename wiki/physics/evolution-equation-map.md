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

Tests:
- `Unit_Tests/unit_test_HLL_flux.c`
- `Unit_Tests/unit_test_hybrid_flux.c`
- `Unit_Tests/unit_test_tabulated_flux.c`
- `Unit_Tests/unit_test_ET_Legacy_HLL_flux.c`
- `Unit_Tests/data_gen/unit_test_data_HLL_flux.c`
- `Unit_Tests/data_gen/unit_test_data_hybrid_flux.c`
- `Unit_Tests/data_gen/unit_test_data_tabulated_flux.c`

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
- NRPy source helpers: `GRHayL/Flux_Source/IGM_All_Source_Terms.py`,
  `GRHayL/Flux_Source/GRHayL_rhs.py`; route generator boundaries through
  [Flux_Source generated NRPy boundary](../gems/flux-source/generated-nrpy-boundary.md)

Tests:
- `Unit_Tests/unit_test_ET_Legacy_flux_source.c`
- `Unit_Tests/data_gen/unit_test_data_ET_Legacy_flux_source.c`
- `Unit_Tests/pert_test_fail_conservatives.c`

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
