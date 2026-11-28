# Evolution Equation Map

> Map from GRMHD conservation laws, stress-energy, fluxes, sources, induction, and Con2Prim recovery in `docs/raw/derivation.md` to C routines and unit tests. · Status: confirmed
> Up: [Physics](index.md)

## Summary

This page maps each GRMHD evolution stage in `docs/raw/derivation.md` and `docs/raw/Induction.dox` to its C implementation files, defining KB leaves, and unit tests without duplicating mathematical derivations.

## Detail

### Conservative Conversion And Stress-Energy Tensor

- **Primitive-to-conservative & limits:** `ghl_compute_conservs` (`compute_conservs.c`), `ghl_compute_conservs_and_Tmunu` (`compute_conservs_and_Tmunu.c`), `ghl_apply_conservative_limits` (`apply_conservative_limits.c`), and `ghl_undensitize_conservatives` (`undensitize_conservatives.c`) implement conservative conversion and limiting ([Limits And Conversions](../gems/con2prim/limits-and-conversions.md)). Tested by `unit_test_compute_conservs_and_Tmunu.c`, `unit_test_apply_conservative_limits.c`, and `unit_test_ET_Legacy_conservs.c`.
- **Stress-energy tensor (`ghl_stress_energy`):** `ghl_compute_smallb_and_b2` (`compute_smallb_and_b2.c`), `ghl_compute_TDNmunu` (`compute_TDNmunu.c`), and `ghl_compute_TUPmunu` (`compute_TUPmunu.c`) compute the magnetic four-vector and combined fluid-EM tensor ([Stress-Energy And `smallb`](../core/stress-energy-smallb.md)). `unit_test_compute_conservs_and_Tmunu.c` exercises `ghl_compute_TDNmunu` and indirectly `smallb` (using comparison helper `pert_test_fail_stress_energy.c` and reference data from `unit_test_data_con2prim_multi_method_hybrid.c`), but does not call `ghl_compute_TUPmunu`.

Evidence (descriptive behavior): `compute_conservs_and_Tmunu.c::ghl_compute_conservs_and_Tmunu`; corroboration: `compute_TDNmunu.c::ghl_compute_TDNmunu`.

### Hydrodynamic Fluxes And Geometric/Neutrino Source Terms

- **Characteristic speeds & HLLE fluxes (`ghl_flux_source.h`):** Direction-split speed kernels (`ghl_calculate_characteristic_speed_dirn{0,1,2}.c`) and direct family/entropy HLLE flux kernels (`hybrid`, `hybrid_entropy`, `tabulated`, `tabulated_entropy`)—both legacy `void` and `_checked` error-returning forms—are the direct public entry points; deprecated generic `ghl_calculate_HLLE_fluxes_dirn*` pointers in `ghl_eos_functions.h` are never assigned by GRHayL ([HLLE Flux Variant Matrix](../gems/flux-source/hlle-flux-variant-matrix.md)). Direct tabulated HLLE definitions remain built without HDF5. Reconstruction (`GRHayL/Reconstruction/`) supplies face states as built public API with no in-repo production caller. Tested by `unit_test_hybrid_flux.c`, `unit_test_tabulated_flux.c`, `unit_test_ET_Legacy_flux_source.c`, and `unit_test_tabulated_eos_compose.c` (note: Induction HLL tests cover magnetic vector-potential fluxes, not hydrodynamic HLLE).
- **Geometric source terms:** `ghl_calculate_source_terms` (`ghl_calculate_source_terms.c`, generated from `generate_flux_source.py`) takes primitives, caller-supplied metric derivatives, and extrinsic curvature, computing tensor algebra inline without accepting `ghl_stress_energy` or calling `compute_TDNmunu`/`compute_TUPmunu` ([Source Terms](../gems/flux-source/source-terms.md)). Tested by `unit_test_ET_Legacy_flux_source.c`.
- **Neutrino matter sources:** `NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c` computes local `R_source` and `Q_source` separately; geometric `ghl_calculate_source_terms` does not add leakage terms ([Physics Model And Sources](../gems/neutrinos/physics-model-and-sources.md)).

Evidence (descriptive behavior): `ghl_calculate_source_terms.c::ghl_calculate_source_terms`; corroboration: `NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c`.

### Induction Vector Potential And Con2Prim Recovery

- **Induction (`ghl_induction.h`):** Separates magnetic HLL fluxes (`HLL_flux_with_B.c`, `HLL_flux_with_Btilde.c`), staggered interpolators (`GRHayL/Induction/Interpolators/`), and Lorenz scalar-potential RHS (`calculate_phitilde_rhs.c`), depending on `Flux_Source` characteristic speeds. Tested by `unit_test_HLL_flux.c`, `unit_test_ET_Legacy_HLL_flux.c`, `unit_test_induction_{ccc_ADM,ccc_BSSN,vvv_ADM}.c`, and `unit_test_ET_Legacy_induction_gauge_rhs.c` (helpers `compute_A_flux_with_B.c`, `compute_A_flux_with_Btilde.c`); see [Induction](../gems/induction/index.md).
- **Con2Prim recovery (`ghl_con2prim.h`):** `con2prim_multi_method.c` dispatches across hybrid (`GRHayL/Con2Prim/Hybrid/`) and tabulated (`GRHayL/Con2Prim/Tabulated/`) solvers with primitive guesses (`guess_primitives.c`), limits (`enforce_primitive_limits_and_compute_u0.c`), and diagnostics (`initialize_diagnostics.c`, where `apply_conservative_limits` accumulates `diagnostics.tau_fix` when raising `tau` to family-specific `tau_atm`). Tested by `unit_test_con2prim_multi_method_hybrid.c`, `unit_test_con2prim_tabulated.c`, `unit_test_con2prim_debug.c`, and `unit_test_hybrid_failure.c`; see [Recovery Flow](../gems/con2prim/recovery-flow.md).

Evidence (descriptive behavior): `calculate_phitilde_rhs.c::ghl_calculate_phitilde_rhs`; corroboration: `apply_conservative_limits.c::ghl_apply_conservative_limits`.

## Sources

- [docs/raw/derivation.md](../../docs/raw/derivation.md) - GRMHD conservative equations, stress-energy, fluxes, and source terms
- [docs/raw/Induction.dox](../../docs/raw/Induction.dox) - vector-potential induction and Lorenz gauge equations
- [GRHayL/Con2Prim/compute_conservs_and_Tmunu.c](../../GRHayL/Con2Prim/compute_conservs_and_Tmunu.c) - combined conservative and `T_{\mu\nu}` calculation
- [GRHayL/Con2Prim/apply_conservative_limits.c](../../GRHayL/Con2Prim/apply_conservative_limits.c) - conservative floor limiting and `tau_fix`
- [GRHayL/GRHayL_Core/compute_TDNmunu.c](../../GRHayL/GRHayL_Core/compute_TDNmunu.c) - covariant stress-energy helper
- [GRHayL/Flux_Source/ghl_calculate_source_terms.c](../../GRHayL/Flux_Source/ghl_calculate_source_terms.c) - inline geometric source kernel
- [GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c](../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c) - local neutrino `R_source`/`Q_source`
- [GRHayL/Induction/calculate_phitilde_rhs.c](../../GRHayL/Induction/calculate_phitilde_rhs.c) - scalar-potential RHS

## See Also

- Depends on: [Variables And Conventions](variables-and-conventions.md)
- Implements: [HLLE Flux Variant Matrix](../gems/flux-source/hlle-flux-variant-matrix.md)
- Implements: [Source Terms](../gems/flux-source/source-terms.md)
- See also: [Recovery Flow](../gems/con2prim/recovery-flow.md)
- See also: [Induction](../gems/induction/index.md)
