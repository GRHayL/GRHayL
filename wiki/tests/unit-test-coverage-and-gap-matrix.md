# Unit Test Coverage And Gap Matrix

> Subsystem-by-subsystem classification of `Unit_Tests/` coverage types, proof limits, and documented test gaps. · Status: confirmed
> Up: [Tests](index.md)

## Summary

This page classifies what GRHayL's `Unit_Tests/` suite exercises directly versus indirectly, where coverage relies on implementation-derived reference binaries or synthetic callbacks, and where explicit test gaps remain. Labels describe repository structure (`direct test`, `replay`, `expected-failure`, `workflow-only`, `manual-only`, `coverage-gap`), never historical run outcomes.

## Detail

### Core, EOS, And Con2Prim Coverage

- **Core (`unit_test_grhayl_core_test_suite.c`, indirect `compute_conservs_and_Tmunu` and `ET_Legacy_conservs`):** Direct checks cover simple/hybrid EOS bounds, simple/hybrid constant atmosphere resets, fixed 4D metric helpers, determinant-enforced ADM init, and valid Con2Prim routine names. **Gaps:** no full `ghl_initialize_*`/`ghl_return_*` round trip, no `ghl_initialize_params` default assertions, no direct clamp (`ghl_imin`/`ghl_imax`/`ghl_iclamp`/`ghl_clamp`), standalone 3D metric helper, or `ghl_compute_TUPmunu` checks, and the tabulated atmosphere block is commented out (`TODO`). See [Core Tests](../core/tests-and-reference-data.md).
- **EOS (`unit_test_piecewise_polytrope.c`, `unit_test_tabulated_eos.c`, `unit_test_tabulated_eos_compose.c`, `unit_test_code_error.c`):** Direct tests cover piecewise polytropes, analytic `simple_table.h5` interpolation/beta-equilibrium, CompOSE Python 100% branch coverage plus C runtime integration, and HDF5/bound error keys. `unit_test_c2p_nn_guess.c` asserts tabulated atmosphere resets indirectly with synthetic metadata and fake callbacks (not table-backed atmosphere qualification). See [EOS Tests](../gems/eos/tests-and-reference-data.md).
- **Con2Prim (`unit_test_con2prim_multi_method_hybrid.c`, `unit_test_con2prim_tabulated.c`, `unit_test_c2p_nn_guess.c`, `unit_test_hybrid_failure.c`, `unit_test_apply_conservative_limits.c`, `unit_test_enforce_primitive_limits_and_compute_u0.c`, `unit_test_compute_conservs_and_Tmunu.c`, `unit_test_con2prim_debug.c`, ET_Legacy):** Covers diagnostics init, conservative/primitive limits, Font1D/Noble1D_entropy2 failure & round-trip behavior, undensitized guesses, and NN HDF5/inference fallback. **Limits:** hybrid reference binaries are implementation-derived regressions; tabulated NN retry remains table/model-dependent; `unit_test_con2prim_debug.c` is manual-only. See [Con2Prim Tests](../gems/con2prim/tests-and-reference-data.md).

Evidence (descriptive behavior): `Unit_Tests/unit_test_grhayl_core_test_suite.c`; corroboration: `Unit_Tests/unit_test_c2p_nn_guess.c`.

### Flux_Source, Induction, Neutrinos, And Reconstruction Coverage

- **Flux_Source (`unit_test_hybrid_flux.c`, `unit_test_tabulated_flux.c`, `unit_test_ET_Legacy_flux_source.c`):** Covers direct hybrid/tabulated HLLE (legacy and `_checked`), callback errors, and ET_Legacy geometric sources. **Limits:** `unit_test_tabulated_flux.c` uses synthetic `test_compute_h_and_cs2.c` (not an EOS-table sound-speed oracle); deprecated generic `ghl_calculate_HLLE_fluxes_dirn*` pointers remain unassigned and untested; `unit_test_HLL_flux.c` belongs to Induction, not Flux_Source. See [Flux_Source Tests](../gems/flux-source/tests-and-reference-data.md).
- **Induction (`unit_test_HLL_flux.c`, `unit_test_ET_Legacy_HLL_flux.c`, `unit_test_induction_{ccc_ADM,ccc_BSSN,vvv_ADM}.c`, `unit_test_ET_Legacy_induction_gauge_rhs.c`):** Covers vector-potential HLL, `ccc_ADM`/`ccc_BSSN`/`vvv_ADM` interpolators, and ET_Legacy gauge RHS. **Gaps:** no standalone `phitilde_rhs`-only unit test and no `vvv_BSSN` interpolation test. See [Induction Tests](../gems/induction/tests-and-reference-data.md).
- **Neutrinos (`unit_test_nrpyleakage_{physics,classifier_fallback,optically_thin_gas,constant_density_sphere,luminosities}.c`, `unit_test_code_error.c`):** Table-free `physics` and `classifier_fallback` verify Fermi branches, nucleon blocking, detailed balance, diffusion, free-nucleon endpoints, chemical-potential reference invariance, asymmetric optical-depth stencil, and portable NaN/Inf classification without HDF5; table-backed replay checks optically thin evolution (with RK lookup failure injection), constant-density sphere, and luminosities on `SLy4`. See [Neutrinos Tests](../gems/neutrinos/tests-and-reference-data.md).
- **Reconstruction (`unit_test_PLM_reconstruction.c`, `unit_test_WENOZ_reconstruction.c`, `unit_test_ET_Legacy_reconstruction.c`):** PLM and WENO-z have direct reference-data tests (note: `.github/run_tests.sh` omits WENO-z while CI workflows run it). **Gap:** no dedicated `unit_test_PPM_reconstruction.c` or PPM data generator exists; PPM is exercised only through `ET_Legacy_reconstruction`. See [Reconstruction Tests](../gems/reconstruction/tests-and-reference-data.md).

Evidence (descriptive behavior): `Unit_Tests/unit_test_nrpyleakage_physics.c`; corroboration: `Unit_Tests/unit_test_tabulated_flux.c`.

### Cross-Cutting Error Keys, Helpers, And Run-Mode Gaps

- **Error keys (`unit_test_code_error.c`):** Expected nonzero exits for keys `0`–`90`; keys `86`–`88` pass enum values directly to `ghl_abort_if_error` without testing producing physics routines; HDF5-only keys skip cleanly under `GHL_DISABLE_HDF5`. See [Expected Failure And Error Keys](expected-failure-and-error-keys.md).
- **Test-only helpers & generators:** Helper C/H files under `Unit_Tests/` and `ghl_unit_tests.h` are not installed production API; `make datagen` compiles generators without running them; ET_Legacy tests do not verify `implementations/GRHayLib/` Cactus integration.

## Sources

- [Unit_Tests/unit_test_grhayl_core_test_suite.c](../../Unit_Tests/unit_test_grhayl_core_test_suite.c) - Core direct suite and commented tabulated block
- [Unit_Tests/unit_test_c2p_nn_guess.c](../../Unit_Tests/unit_test_c2p_nn_guess.c) - NN inference/HDF5 and synthetic atmosphere fallback checks
- [Unit_Tests/unit_test_tabulated_flux.c](../../Unit_Tests/unit_test_tabulated_flux.c) - tabulated HLLE test using synthetic callback
- [Unit_Tests/unit_test_nrpyleakage_physics.c](../../Unit_Tests/unit_test_nrpyleakage_physics.c) - table-free NRPyLeakage physics checks
- [Unit_Tests/unit_test_code_error.c](../../Unit_Tests/unit_test_code_error.c) - expected-failure error-key suite
- [.github/run_tests.sh](../../.github/run_tests.sh) - broad runner coverage vs. CI matrices

## See Also

- Depends on: [Unit Test Catalog](unit-test-catalog.md)
- Depends on: [Data Generators And Helpers](data-generators-and-helpers.md)
- See also: [Runner And Build Outputs](runner-and-build-outputs.md)
- See also: [ET Legacy Comparison](et-legacy-comparison.md)
