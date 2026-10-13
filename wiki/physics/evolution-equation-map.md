# Evolution Equation Map

Use [docs/raw/derivation.md](../../docs/raw/derivation.md) for equation authority. This page maps equation concepts to code and tests without repeating the derivation.

## Conservative Equations

Concept:
- Evolution is written for conservative variables rather than primitive variables.
- The derivation summarizes density, energy, momentum, and magnetic-field conservative equations.

Code map:
- Conservative struct: [GRHayL/include/ghl.h](../../GRHayL/include/ghl.h)
- Primitive-to-conservative conversion: [GRHayL/Con2Prim/compute_conservs.c](../../GRHayL/Con2Prim/compute_conservs.c)
- Conversion plus stress-energy: [GRHayL/Con2Prim/compute_conservs_and_Tmunu.c](../../GRHayL/Con2Prim/compute_conservs_and_Tmunu.c)
- Conservative limiting: [GRHayL/Con2Prim/apply_conservative_limits.c](../../GRHayL/Con2Prim/apply_conservative_limits.c)
- Undensitization for recovery: [GRHayL/Con2Prim/undensitize_conservatives.c](../../GRHayL/Con2Prim/undensitize_conservatives.c)
- Con2Prim helper route: [limits and conversions](../gems/con2prim/limits-and-conversions.md)

Tests:
- [Unit_Tests/unit_test_compute_conservs_and_Tmunu.c](../../Unit_Tests/unit_test_compute_conservs_and_Tmunu.c)
- [Unit_Tests/unit_test_apply_conservative_limits.c](../../Unit_Tests/unit_test_apply_conservative_limits.c)
- [Unit_Tests/unit_test_ET_Legacy_conservs.c](../../Unit_Tests/unit_test_ET_Legacy_conservs.c)
- [Unit_Tests/data_gen/unit_test_data_ET_Legacy_conservs.c](../../Unit_Tests/data_gen/unit_test_data_ET_Legacy_conservs.c)

## Stress-Energy Tensor

Concept:
- `docs/raw/derivation.md` derives the electromagnetic contribution and combined stress-energy tensor with GRHayL magnetic rescaling.
- Stress-energy feeds conservative energy, momentum, fluxes, and source terms.

Code map:
- Public struct: `ghl_stress_energy` in `GRHayL/include/ghl.h`
- Small-b and magnetic scalar: [GRHayL/GRHayL_Core/compute_smallb_and_b2.c](../../GRHayL/GRHayL_Core/compute_smallb_and_b2.c)
- Tensor helpers: [GRHayL/GRHayL_Core/compute_TDNmunu.c](../../GRHayL/GRHayL_Core/compute_TDNmunu.c), [GRHayL/GRHayL_Core/compute_TUPmunu.c](../../GRHayL/GRHayL_Core/compute_TUPmunu.c)
- Conservative and tensor combined path: `GRHayL/Con2Prim/compute_conservs_and_Tmunu.c`

Tests:
- [Unit_Tests/pert_test_fail_stress_energy.c](../../Unit_Tests/pert_test_fail_stress_energy.c)
- `Unit_Tests/unit_test_compute_conservs_and_Tmunu.c`
- Fixture producer:
  [Unit_Tests/data_gen/unit_test_data_con2prim_multi_method_hybrid.c](../../Unit_Tests/data_gen/unit_test_data_con2prim_multi_method_hybrid.c).

`pert_test_fail_stress_energy.c` is a comparison helper, not a standalone
test. The combined test also exercises `ghl_compute_TDNmunu`; smallb is
covered indirectly, while it does not call `ghl_compute_TUPmunu`. Follow the
[Core stress-energy contract](../core/stress-energy-smallb-contract.md) for
component-write and coverage limits.

## Flux Terms

Concept:
- Flux terms move conservative quantities through cell faces.
- HLLE flux paths are split by direction, EOS family, and entropy support.
- Reconstruction provides face values before flux evaluation.

Code map:
- Public flux API: [GRHayL/include/ghl_flux_source.h](../../GRHayL/include/ghl_flux_source.h)
- Characteristic speeds: [GRHayL/Flux_Source/ghl_calculate_characteristic_speed_dirn0.c](../../GRHayL/Flux_Source/ghl_calculate_characteristic_speed_dirn0.c), [GRHayL/Flux_Source/ghl_calculate_characteristic_speed_dirn1.c](../../GRHayL/Flux_Source/ghl_calculate_characteristic_speed_dirn1.c), [GRHayL/Flux_Source/ghl_calculate_characteristic_speed_dirn2.c](../../GRHayL/Flux_Source/ghl_calculate_characteristic_speed_dirn2.c); route caller contract through [Flux_Source characteristic speeds](../gems/flux-source/characteristic-speeds-contract.md)
- Hybrid, hybrid entropy, tabulated, and tabulated entropy HLLE variants: route
  direction/EOS variant mapping through [Flux_Source HLLE flux variant matrix](../gems/flux-source/hlle-flux-variant-matrix.md)
- Reconstruction sources: [GRHayL/Reconstruction/](../../GRHayL/Reconstruction)

Availability boundary:

- Family-specific direct HLLE functions, legacy `void` names and their
  `_checked` counterparts, are the direct public HLLE entry points. Deprecated
  generic `ghl_calculate_HLLE_fluxes_dirn0/1/2` pointers remain declared in
  [GRHayL/include/ghl_eos_functions.h](../../GRHayL/include/ghl_eos_functions.h) and exported, but EOS initialization
  never assigns them. A caller can supply an exact-signature compatibility
  callback; direct routines have different signatures and cannot simply be
  assigned to those pointers. See the
  [HLLE variant matrix](../gems/flux-source/hlle-flux-variant-matrix.md).
- No-HDF5 builds retain the direct tabulated HLLE definitions; their
  table-dependent tests and generators remain excluded.
- Reconstruction is built public caller API, but no production source caller is
  visible in this repository. Tests demonstrate family kernels with documented
  assertion/generator limits; they do not prove an evolution integration.

Tests:
- [Unit_Tests/unit_test_hybrid_flux.c](../../Unit_Tests/unit_test_hybrid_flux.c)
- [Unit_Tests/unit_test_tabulated_flux.c](../../Unit_Tests/unit_test_tabulated_flux.c)
- [Unit_Tests/unit_test_ET_Legacy_flux_source.c](../../Unit_Tests/unit_test_ET_Legacy_flux_source.c)
- [Unit_Tests/unit_test_tabulated_eos_compose.c](../../Unit_Tests/unit_test_tabulated_eos_compose.c)
- [Unit_Tests/data_gen/unit_test_data_hybrid_flux.c](../../Unit_Tests/data_gen/unit_test_data_hybrid_flux.c)
- [Unit_Tests/data_gen/unit_test_data_tabulated_flux.c](../../Unit_Tests/data_gen/unit_test_data_tabulated_flux.c)

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
- Source implementation: [GRHayL/Flux_Source/ghl_calculate_source_terms.c](../../GRHayL/Flux_Source/ghl_calculate_source_terms.c);
  route caller-owned metric derivative and extrinsic-curvature contract through
  [Flux_Source source-term contract](../gems/flux-source/source-terms-contract.md)
- Related tensor helpers: `GRHayL/GRHayL_Core/compute_TDNmunu.c`, `GRHayL/GRHayL_Core/compute_TUPmunu.c`.
  The source kernel takes primitives and computes its tensor algebra inline;
  it accepts no precomputed stress-energy struct and calls neither helper.
- NRPy 2 source generator: [GRHayL/Flux_Source/generate_flux_source.py](../../GRHayL/Flux_Source/generate_flux_source.py);
  route generation through [GRHayL/Flux_Source/generate_flux_source.sh](../../GRHayL/Flux_Source/generate_flux_source.sh) and
  [Flux_Source generated NRPy boundary](../gems/flux-source/generated-nrpy-boundary.md)

Tests:
- `Unit_Tests/unit_test_ET_Legacy_flux_source.c`
- [Unit_Tests/data_gen/unit_test_data_ET_Legacy_flux_source.c](../../Unit_Tests/data_gen/unit_test_data_ET_Legacy_flux_source.c)
- [Unit_Tests/pert_test_fail_conservatives.c](../../Unit_Tests/pert_test_fail_conservatives.c)

Neutrino matter-source contributions are a separate route:
[NRPyLeakage source terms](../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c)
return local `R_source` and `Q_source`. Follow the
[Neutrinos physics contract](../gems/neutrinos/physics-and-eos-contract.md)
for units, signs, and the caller-owned conservative/metric conversion. The
geometric Flux_Source kernel does not incorporate those leakage outputs.

## Induction and Vector Potential

Concept:
- Magnetic evolution uses vector potential and densitized scalar potential paths documented in [docs/raw/Induction.dox](../../docs/raw/Induction.dox).
- The induction gem separates magnetic HLL flux terms from gauge/scalar-potential RHS terms.
- Keep equation authority in `docs/raw/Induction.dox` and `docs/raw/derivation.md`; use KB pages only for routing.

Code map:
- Public induction API: [GRHayL/include/ghl_induction.h](../../GRHayL/include/ghl_induction.h)
- Magnetic HLL flux: [GRHayL/Induction/HLL_flux_with_B.c](../../GRHayL/Induction/HLL_flux_with_B.c), [GRHayL/Induction/HLL_flux_with_Btilde.c](../../GRHayL/Induction/HLL_flux_with_Btilde.c); route caller contract through [Induction HLL flux contract](../gems/induction/hll-flux-contract.md)
- Interpolation: [GRHayL/Induction/Interpolators/](../../GRHayL/Induction/Interpolators)
- Scalar-potential RHS: [GRHayL/Induction/calculate_phitilde_rhs.c](../../GRHayL/Induction/calculate_phitilde_rhs.c); route gauge/scalar RHS contract through [Induction gauge RHS contract](../gems/induction/gauge-rhs-contract.md)
- Characteristic speed dependency:
  [Flux_Source characteristic speeds](../gems/flux-source/characteristic-speeds-contract.md)

Tests:
- [Unit_Tests/unit_test_induction_ccc_ADM.c](../../Unit_Tests/unit_test_induction_ccc_ADM.c)
- [Unit_Tests/unit_test_induction_ccc_BSSN.c](../../Unit_Tests/unit_test_induction_ccc_BSSN.c)
- [Unit_Tests/unit_test_induction_vvv_ADM.c](../../Unit_Tests/unit_test_induction_vvv_ADM.c)
- [Unit_Tests/unit_test_ET_Legacy_induction_gauge_rhs.c](../../Unit_Tests/unit_test_ET_Legacy_induction_gauge_rhs.c)
- [Unit_Tests/unit_test_HLL_flux.c](../../Unit_Tests/unit_test_HLL_flux.c)
- [Unit_Tests/unit_test_ET_Legacy_HLL_flux.c](../../Unit_Tests/unit_test_ET_Legacy_HLL_flux.c)
- [Unit_Tests/data_gen/unit_test_data_HLL_flux.c](../../Unit_Tests/data_gen/unit_test_data_HLL_flux.c)
- [Unit_Tests/compute_A_flux_with_B.c](../../Unit_Tests/compute_A_flux_with_B.c)
- [Unit_Tests/compute_A_flux_with_Btilde.c](../../Unit_Tests/compute_A_flux_with_Btilde.c)
- [Unit_Tests/data_gen/unit_test_data_induction_interpolation.c](../../Unit_Tests/data_gen/unit_test_data_induction_interpolation.c)
- Fixture and run routes: [Induction tests and fixtures](../gems/induction/tests-and-fixtures.md),
  [Induction verification workflows](../gems/induction/verification-workflows.md)

## Conservative-to-Primitive Recovery

Concept:
- Conservative evolution requires recovering primitive variables after updates.
- Recovery is numerical and EOS-dependent.

Code map:
- Public recovery API: [GRHayL/include/ghl_con2prim.h](../../GRHayL/include/ghl_con2prim.h)
- Multi-method dispatch: [GRHayL/Con2Prim/con2prim_multi_method.c](../../GRHayL/Con2Prim/con2prim_multi_method.c)
- Hybrid methods: [GRHayL/Con2Prim/Hybrid/](../../GRHayL/Con2Prim/Hybrid)
- Tabulated methods: [GRHayL/Con2Prim/Tabulated/](../../GRHayL/Con2Prim/Tabulated)
- Guess and bounds: [GRHayL/Con2Prim/guess_primitives.c](../../GRHayL/Con2Prim/guess_primitives.c), [GRHayL/Con2Prim/enforce_primitive_limits_and_compute_u0.c](../../GRHayL/Con2Prim/enforce_primitive_limits_and_compute_u0.c)
- Diagnostics: [GRHayL/Con2Prim/initialize_diagnostics.c](../../GRHayL/Con2Prim/initialize_diagnostics.c)
- KB routes: [solver matrix](../gems/con2prim/solver-matrix.md),
  [recovery flow](../gems/con2prim/recovery-flow.md),
  [limits and conversions](../gems/con2prim/limits-and-conversions.md)

Tests:
- [Unit_Tests/unit_test_con2prim_multi_method_hybrid.c](../../Unit_Tests/unit_test_con2prim_multi_method_hybrid.c)
- [Unit_Tests/unit_test_con2prim_tabulated.c](../../Unit_Tests/unit_test_con2prim_tabulated.c)
- [Unit_Tests/unit_test_con2prim_debug.c](../../Unit_Tests/unit_test_con2prim_debug.c)
- [Unit_Tests/unit_test_hybrid_failure.c](../../Unit_Tests/unit_test_hybrid_failure.c)
- [Unit_Tests/data_gen/unit_test_data_con2prim_multi_method_hybrid.c](../../Unit_Tests/data_gen/unit_test_data_con2prim_multi_method_hybrid.c)
- Fixture route: [Con2Prim tests and fixtures](../gems/con2prim/tests-and-fixtures.md)

Recovery floor note: `apply_conservative_limits` sets and accumulates
`diagnostics->tau_fix` when raising tau to `tau_atm`; family-dependent
`tau_atm` setup is routed through [EOS dispatch](../core/eos-dispatch-contract.md).
