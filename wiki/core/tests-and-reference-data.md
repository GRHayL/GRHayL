# Core Tests And Reference Data

> Direct Core unit tests, reference-data generators, shared test helpers, indirect Core-dependent suites, and Core coverage gaps. · Status: confirmed
> Up: [GRHayL Core](index.md)

## Summary

`Unit_Tests/unit_test_grhayl_core_test_suite.c` directly tests simple/hybrid EOS bounds, simple/hybrid constant atmosphere resets, fixed 4D metric helpers, determinant-enforced ADM initialization (via `grhayl_core_test_suite_input.bin`), and valid `ghl_get_con2prim_routine_name` keys, while stress-energy, `u0` limiting, and error keys are covered by adjacent suites. Several Core helpers (`ghl_compute_TUPmunu`, clamp utilities, standalone 3D metric helpers, `ghl_initialize_params` defaults, and tabulated atmosphere setup) lack direct Core-suite assertions.

## Detail

### Direct Core Tests And Reference-Data Generation

- **`unit_test_grhayl_core_test_suite.c`:** checks simple EOS rectangular extrema and defaults, simple/hybrid zero-density floors, hybrid EOS setup, `ghl_set_prims_to_constant_atm` for simple/hybrid (rejecting nonfinite outputs before tolerances), fixed-array `ghl_compute_vec2_from_vec4D` and `ghl_raise_lower_vector_4D`, determinant-enforced ADM initialization from `grhayl_core_test_suite_input.bin`, and valid `ghl_get_con2prim_routine_name` keys. Its tabulated EOS atmosphere block is commented out (`TODO`), and its `Y_e`/temperature atmosphere assertions sit behind an unreached `eos_type == 2` guard.
- **`data_gen/unit_test_data_grhayl_core_test_suite.c`:** calls `ghl_randomize_metric` to write lowercase `grhayl_core_test_suite_input.bin`.
- **`unit_test_code_error.c`:** covers invalid EOS bounds, singular `u0` (key `4`), invalid `ghl_get_con2prim_routine_name` (key `33`), and HDF5/table/bracketing errors as expected nonzero process exits.

Evidence (descriptive behavior): `unit_test_grhayl_core_test_suite.c::grhayl_core_test_suite_input.bin`; corroboration: `unit_test_data_grhayl_core_test_suite.c`.

### Indirect Core-Dependent Suites, Helpers, And Coverage Gaps

- **Stress-energy & conservatives:** `unit_test_compute_conservs_and_Tmunu.c` (with `pert_test_fail_stress_energy.c` and Con2Prim binary data `metric_Bfield_initial_data.bin`, `compute_conservs_and_Tmunu_*.bin`) checks combined conservative/stress-energy output and standalone `ghl_compute_TDNmunu`, exercising `ghl_compute_smallb_and_b2` indirectly without testing `ghl_compute_TUPmunu`.
- **Velocity & legacy repacking:** `unit_test_enforce_primitive_limits_and_compute_u0.c` tests Con2Prim limiting plus focused direct `ghl_limit_v_and_compute_u0` cap cases; `unit_test_ET_Legacy_conservs.c` indirectly exercises metric/primitive/conservative/stress-energy packing and `ghl_return_*` repacking; Con2Prim and flux tests repeatedly call Core initializers (flux face states stay below the Lorentz cap).
- **Helper-only files:** `ghl_unit_tests.h` (not installed), `pert_test_fail_{conservatives,primitives,stress_energy}.c`, and `randomize_{metric,primitives}.c` (`randomize_metric.c` fixes `lapse = 1.0`, `shift = 0`); see [Data Generators And Helpers](../tests/data-generators-and-helpers.md).
- **Direct coverage gaps:** no complete `ghl_initialize_*`/`ghl_return_*` round-trip test; no assertions on `ghl_initialize_params` Con2Prim/PPM defaults; no direct test of uninitialized EOS pointers, full dispatch-target tables, or invalid `ghl_eos_t`; no Core-suite check of singular simple `Gamma == 1`, clamp functions (`ghl_imin`/`ghl_imax`/`ghl_iclamp`/`ghl_clamp`), standalone 3D metric helpers, `ghl_compute_TUPmunu`, or tabulated atmosphere resets.

Evidence (descriptive behavior): `unit_test_compute_conservs_and_Tmunu.c::ghl_compute_TDNmunu`; corroboration: `randomize_metric.c`.

## Sources

- [Unit_Tests/unit_test_grhayl_core_test_suite.c](../../Unit_Tests/unit_test_grhayl_core_test_suite.c) - direct Core test suite
- [Unit_Tests/data_gen/unit_test_data_grhayl_core_test_suite.c](../../Unit_Tests/data_gen/unit_test_data_grhayl_core_test_suite.c) - Core metric reference generator
- [Unit_Tests/unit_test_code_error.c](../../Unit_Tests/unit_test_code_error.c) - expected-failure error-key suite
- [Unit_Tests/unit_test_compute_conservs_and_Tmunu.c](../../Unit_Tests/unit_test_compute_conservs_and_Tmunu.c) - `ghl_compute_TDNmunu` and conservative checks
- [Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c](../../Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c) - direct `ghl_limit_v_and_compute_u0` cases
- [Unit_Tests/unit_test_ET_Legacy_conservs.c](../../Unit_Tests/unit_test_ET_Legacy_conservs.c) - indirect packing/repacking coverage
- [Unit_Tests/randomize_metric.c](../../Unit_Tests/randomize_metric.c) - random metric helper
- [Unit_Tests/pert_test_fail_stress_energy.c](../../Unit_Tests/pert_test_fail_stress_energy.c) - stress-energy comparator

## See Also

- Depends on: [Data Generators And Helpers](../tests/data-generators-and-helpers.md)
- See also: [Unit Test Catalog](../tests/unit-test-catalog.md)
- See also: [Con2Prim Tests And Reference Data](../gems/con2prim/tests-and-reference-data.md)
- See also: [EOS Tests And Reference Data](../gems/eos/tests-and-reference-data.md)
