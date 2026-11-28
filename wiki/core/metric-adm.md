# Metric And ADM Helpers

> Core ADM 3-metric initialization, determinant enforcement, 4-metric auxiliary construction, and inline `ghl_metric_helpers.h` vector operations. · Status: contested
> Up: [GRHayL Core](index.md)

## Summary

Core provides `ghl_initialize_metric`, `ghl_compute_ADM_auxiliaries`, `ghl_enforce_detgtij_and_initialize_ADM_metric`, and header-only 3D/4D vector raise/lower and inner-product helpers in `ghl_metric_helpers.h`. `ghl_compute_vec2_from_vec3D` and `ghl_compute_vec2_from_vec4D` return the `double` contraction; the others return `void`. None validates inputs, and current determinant handling divides by `fabs(detgamma)` while preserving negative determinant signs after warning.

## Detail

### Metric Initialization, ADM Auxiliaries, And Determinant Enforcement

Claim status: contested; contradiction: CONTR-0011. See [CONTR-0011](../contradictions.md#contr-0011).

- **`ghl_initialize_metric` (`initialize_metric.c`):** copies lapse, shift, and six symmetric 3-metric components into `ghl_metric_quantities`, then computes `lapseinv`, `lapseinv2`, `detgamma`, `sqrt_detgamma`, and `gammaUU`. It stores `fabs(raw_determinant)` before computing `sqrt_detgamma` and divides cofactors by `fabs(raw_determinant)`—giving a true inverse only for positive-determinant spatial metrics; zero, negative, or nonfinite inputs are unchecked.
- **`ghl_compute_ADM_auxiliaries` (`compute_ADM_auxiliaries.c`):** lowers shift with `ghl_raise_lower_vector_3D` and populates 4-metric `g4DD` and inverse `g4UU` in `ghl_ADM_aux_quantities`, assuming `metric_adm` already holds a consistent metric and inverse.
- **`ghl_enforce_detgtij_and_initialize_ADM_metric` (`enforce_detgij_and_initialize_ADM_metric.c`):** calls `ghl_initialize_metric`, computes a local `ghl_compute_ADM_auxiliaries` result that it never reads or returns, forms conformal metric components, scales by `fabs(1.0 / cbrt(gtijdet))`, rebuilds the ADM 3-metric, and re-invokes `ghl_initialize_metric`. For negative `gtijdet`, `fabs(...)` preserves the negative sign rather than restoring `+1`, emits a warning, and returns `void` with no error status.

Evidence (descriptive behavior): `initialize_metric.c::ghl_initialize_metric`; corroboration: `enforce_detgij_and_initialize_ADM_metric.c::ghl_enforce_detgtij_and_initialize_ADM_metric`.

### Header-Only Inline Vector Helpers And Test Coverage

`GRHayL/include/ghl_metric_helpers.h` defines inline routines assuming symmetric metric arrays (reading only the diagonal plus one triangle) with no null, symmetry, or finiteness checks:
- `ghl_raise_lower_vector_4D` and `ghl_raise_lower_vector_3D`: output arrays must **not** alias input vectors because output writes occur while later components still read the input.
- `ghl_compute_vec2_from_vec4D` and `ghl_compute_vec2_from_vec3D`: compute squared vector norms.

`Unit_Tests/unit_test_grhayl_core_test_suite.c` (using `grhayl_core_test_suite_input.bin` from `unit_test_data_grhayl_core_test_suite.c`) directly tests `ghl_compute_vec2_from_vec4D`, `ghl_raise_lower_vector_4D`, and `ghl_enforce_detgtij_and_initialize_ADM_metric` (including a negative-determinant warning trigger without asserting positive conformal determinant). Because the determinant-enforced path discards its local ADM auxiliaries, `ghl_compute_ADM_auxiliaries` is numerically checked indirectly via `unit_test_compute_conservs_and_Tmunu.c`, and no standalone fixed-array check for `ghl_compute_vec2_from_vec3D` exists in the Core suite.

Evidence (descriptive behavior): `ghl_metric_helpers.h::ghl_raise_lower_vector_4D`; corroboration: `unit_test_grhayl_core_test_suite.c`.

## Sources

- [GRHayL/include/ghl.h](../../GRHayL/include/ghl.h) - `ghl_metric_quantities`, `ghl_ADM_aux_quantities`
- [GRHayL/include/ghl_metric_helpers.h](../../GRHayL/include/ghl_metric_helpers.h) - 3D/4D raise/lower and vector-square inline helpers
- [GRHayL/GRHayL_Core/initialize_metric.c](../../GRHayL/GRHayL_Core/initialize_metric.c) - `ghl_initialize_metric`
- [GRHayL/GRHayL_Core/compute_ADM_auxiliaries.c](../../GRHayL/GRHayL_Core/compute_ADM_auxiliaries.c) - `ghl_compute_ADM_auxiliaries`
- [GRHayL/GRHayL_Core/enforce_detgij_and_initialize_ADM_metric.c](../../GRHayL/GRHayL_Core/enforce_detgij_and_initialize_ADM_metric.c) - `ghl_enforce_detgtij_and_initialize_ADM_metric`
- [Unit_Tests/unit_test_grhayl_core_test_suite.c](../../Unit_Tests/unit_test_grhayl_core_test_suite.c) - direct 4D helper and determinant-enforced tests
- [Unit_Tests/data_gen/unit_test_data_grhayl_core_test_suite.c](../../Unit_Tests/data_gen/unit_test_data_grhayl_core_test_suite.c) - Core metric reference generator
- [Unit_Tests/unit_test_compute_conservs_and_Tmunu.c](../../Unit_Tests/unit_test_compute_conservs_and_Tmunu.c) - indirect ADM-auxiliary numerical checks

## See Also

- Depends on: [Shared Parameters And Enums](shared-parameters-and-enums.md)
- Validated by: [Core Tests And Reference Data](tests-and-reference-data.md)
- See also: [Contradictions](../contradictions.md#contr-0011)
- See also: [Variables And Conventions](../physics/variables-and-conventions.md)
