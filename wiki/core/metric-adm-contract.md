# Metric/ADM Contract

This page routes Core metric packing, ADM auxiliary construction, and metric
inline-helper questions. Ground truth is `GRHayL/include/ghl.h`,
`GRHayL/include/ghl_metric_helpers.h`, the Core metric source files named
below, and the Core test files named below.

## Metric Packing

- `ghl_initialize_metric` copies lapse, shift, and six symmetric 3-metric
  inputs into `ghl_metric_quantities`, then computes `lapseinv`,
  `lapseinv2`, `detgamma`, `sqrt_detgamma`, and `gammaUU`
  (`GRHayL/GRHayL_Core/initialize_metric.c`).
- The determinant stored by `ghl_initialize_metric` is absolute-valued before
  `sqrt_detgamma` is computed (`GRHayL/GRHayL_Core/initialize_metric.c`).
- Public declarations and the `ghl_metric_quantities` layout live in
  `GRHayL/include/ghl.h`; the Doxygen group route is
  `docs/raw/GRHayL_Core.dox`.

This initializer returns `void` and does not validate inputs. Its direct
division requires nonzero lapse and nonzero determinant. The cofactor inverse
is divided by `fabs(raw_determinant)`, not the signed determinant, so it is an
actual inverse only for the expected positive-determinant spatial metric.
Zero, negative, non-finite, or non-positive-definite metric inputs are not a
checked failure mode (`GRHayL/GRHayL_Core/initialize_metric.c`).

## ADM Auxiliaries

- `ghl_compute_ADM_auxiliaries` builds the covariant 4-metric `g4DD` and
  inverse 4-metric `g4UU` in `ghl_ADM_aux_quantities`
  (`GRHayL/GRHayL_Core/compute_ADM_auxiliaries.c`).
- It lowers the shift through `ghl_raise_lower_vector_3D`, then uses the ADM
  lapse, inverse lapse, shift, `gammaDD`, and `gammaUU` to populate the 4D
  metric arrays (`GRHayL/GRHayL_Core/compute_ADM_auxiliaries.c`;
  `GRHayL/include/ghl_metric_helpers.h`).
- It returns `void` and assumes `metric_adm` already contains a mutually
  consistent ADM metric and inverse. It performs no null, finiteness,
  invertibility, or signature checks
  (`GRHayL/GRHayL_Core/compute_ADM_auxiliaries.c`).

## Determinant-Enforced Path

- `ghl_enforce_detgtij_and_initialize_ADM_metric` first calls
  `ghl_initialize_metric`, then `ghl_compute_ADM_auxiliaries`, computes BSSN
  conformal metric components, normalizes the magnitude of the conformal
  determinant, rebuilds the ADM 3-metric, and calls `ghl_initialize_metric`
  again
  (`GRHayL/GRHayL_Core/enforce_detgij_and_initialize_ADM_metric.c`).
- The ADM-auxiliary result in this routine is local and is not returned; the
  current source computes it but does not subsequently read it.
- For a negative `gtijdet`, the `fabs(1/cbrt(gtijdet))` factor preserves the
  negative sign rather than producing `detgtij = +1`. The routine rebuilds the
  metric, emits a warning, returns `void`, and supplies no error status
  (`GRHayL/GRHayL_Core/enforce_detgij_and_initialize_ADM_metric.c`).

## Inline Helpers

`GRHayL/include/ghl_metric_helpers.h` defines four header-only helpers:

- `ghl_raise_lower_vector_4D` raises or lowers a 4-vector with caller-provided
  4D metric/inverse metric arrays.
- `ghl_raise_lower_vector_3D` raises or lowers a 3-vector with caller-provided
  3D metric/inverse metric arrays.
- `ghl_compute_vec2_from_vec4D` computes a 4-vector square from a 4D
  metric/inverse metric array.
- `ghl_compute_vec2_from_vec3D` computes a 3-vector square from a 3D
  metric/inverse metric array.

All four helpers assume symmetric metric arrays and read the diagonal plus one
triangle; they do not verify symmetry. Raise/lower outputs must not alias their
input vector because assignments occur while later expressions still read the
input. No helper checks pointers, dimensions, finiteness, or metric validity
(`GRHayL/include/ghl_metric_helpers.h`).

## Test Coverage Routes

- `Unit_Tests/unit_test_grhayl_core_test_suite.c` directly checks
  `ghl_compute_vec2_from_vec4D` and `ghl_raise_lower_vector_4D` with fixed
  arrays.
- The same Core suite checks
  `ghl_enforce_detgtij_and_initialize_ADM_metric` against generated metric
  fixtures. It compares returned metric fields with inputs, including after it
  perturbs a metric to trigger the negative-determinant warning; it does not
  assert a positive conformal determinant or a failure code
  (`Unit_Tests/unit_test_grhayl_core_test_suite.c`).
- `ghl_raise_lower_vector_3D` is indirectly covered through
  `ghl_compute_ADM_auxiliaries` because that routine calls the 3D helper before
  building 4D metric arrays (`GRHayL/GRHayL_Core/compute_ADM_auxiliaries.c`;
  `Unit_Tests/unit_test_grhayl_core_test_suite.c`).
- No direct fixed-array test for `ghl_compute_vec2_from_vec3D` appears in the
  Core suite; route new direct 3D-helper tests to
  `Unit_Tests/unit_test_grhayl_core_test_suite.c`.
- Fixture naming drift is source/test-confirmed: the generator opens
  `grhayL_core_test_suite_input.bin`, while the suite reads
  `grhayl_core_test_suite_input.bin`
  (`Unit_Tests/data_gen/unit_test_data_grhayl_core_test_suite.c`;
  `Unit_Tests/unit_test_grhayl_core_test_suite.c`).
