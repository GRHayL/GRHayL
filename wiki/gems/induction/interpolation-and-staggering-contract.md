# Induction Interpolation And Staggering Contract

## Routing Purpose

Use this page when checking Induction interpolator API shape, vector-potential
staggering assumptions, output fields, or test coverage. Start from the Induction
hub [`wiki/gems/induction.md`](../induction.md), then verify public signatures in
[`GRHayL/include/ghl_induction.h`](../../../GRHayL/include/ghl_induction.h) and
implementation behavior in
[`GRHayL/Induction/Interpolators/`](../../../GRHayL/Induction/Interpolators/).

Repo-local source remains authority. GRHayL receives caller-packed stencils, not
grid metadata, so callers own correct centering, indexing, and ghost-zone reach
before calling these wrappers.

## Public API

[`ghl_induction.h`](../../../GRHayL/include/ghl_induction.h) declares three
public induction interpolators:

- `ghl_interpolate_with_cell_centered_ADM`: source
  [`interpolate_with_cell_centered_ADM.c`](../../../GRHayL/Induction/Interpolators/interpolate_with_cell_centered_ADM.c),
  test helper [`compute_ccc_ADM.c`](../../../Unit_Tests/compute_ccc_ADM.c), and
  test [`unit_test_induction_ccc_ADM.c`](../../../Unit_Tests/unit_test_induction_ccc_ADM.c).
- `ghl_interpolate_with_cell_centered_BSSN`: source
  [`interpolate_with_cell_centered_BSSN.c`](../../../GRHayL/Induction/Interpolators/interpolate_with_cell_centered_BSSN.c),
  test helper [`compute_ccc_BSSN.c`](../../../Unit_Tests/compute_ccc_BSSN.c), and
  test [`unit_test_induction_ccc_BSSN.c`](../../../Unit_Tests/unit_test_induction_ccc_BSSN.c).
- `ghl_interpolate_with_vertex_centered_ADM`: source
  [`interpolate_with_vertex_centered_ADM.c`](../../../GRHayL/Induction/Interpolators/interpolate_with_vertex_centered_ADM.c),
  test helper [`compute_vvv_ADM.c`](../../../Unit_Tests/compute_vvv_ADM.c), and
  test [`unit_test_induction_vvv_ADM.c`](../../../Unit_Tests/unit_test_induction_vvv_ADM.c).

The public output type is `ghl_induction_interp_vars`. Its fields are:

- `alpha`
- `alpha_Phi_minus_betaj_A_j`
- `sqrtg_Ai[3]`
- `betai[3]`

The cell-centered ADM and BSSN wrappers fill all of these fields for later
scalar-potential RHS work and for the gauge contribution to the
vector-potential RHS. The vertex-centered ADM wrapper fills `sqrtg_Ai[3]` and
`alpha_Phi_minus_betaj_A_j` only; callers must not read `alpha` or `betai[3]`
from that wrapper. The wrappers do not infer grid locations from metadata.

## Public Stencil Shapes

The public signatures encode these stencil shapes:

| Input | Shape | Used by |
| --- | --- | --- |
| `metric_stencil` or `metric` | `2x2x2` `ghl_metric_quantities` | all three wrappers |
| `psi_stencil` | `2x2x2` scalar `psi` | `ghl_interpolate_with_cell_centered_BSSN` only |
| `Ax_stencil` | `3x3x3` scalar `A_x` | all three wrappers |
| `Ay_stencil` | `3x3x3` scalar `A_y` | all three wrappers |
| `Az_stencil` | `3x3x3` scalar `A_z` | all three wrappers |
| `phitilde` | scalar `tildePhi` value | all three wrappers |

The source comments describe exact cell offsets for these arrays. This KB page
keeps only route-level contract: `A_i` is staggered in directions perpendicular
to its component, and `tildePhi` is fully staggered.

## Public Wrappers And Internal Helpers

The public wrappers live in `ghl_induction.h` and are compiled through
[`GRHayL/Induction/Interpolators/make.code.defn`](../../../GRHayL/Induction/Interpolators/make.code.defn).
They call helper routines declared in
[`ghl_induction_helpers.h`](../../../GRHayL/Induction/Interpolators/ghl_induction_helpers.h)
and implemented in
[`interpolate_helper.c`](../../../GRHayL/Induction/Interpolators/interpolate_helper.c):

- `ghl_A_i_avg` averages `A_i` components to `tildePhi`, `A_x`, `A_y`, and
  `A_z` locations.
- `ghl_ADM_cell_interp` handles cell-centered ADM metric interpolation.
- `ghl_BSSN_cell_interp` handles cell-centered BSSN metric and `psi_stencil`
  interpolation.
- `ghl_ADM_vertex_interp` handles vertex-centered ADM metric interpolation.

Treat these helpers as internal interpolation machinery. Callers should use the
public wrappers unless they are deliberately working inside the Induction
interpolator implementation.

## Coverage

Interpolation coverage is fixture-based:

- `ccc_ADM` calls `ghl_interpolate_with_cell_centered_ADM` through
  [`compute_ccc_ADM.c`](../../../Unit_Tests/compute_ccc_ADM.c). The unit test
  compares `alpha`, all shift components, `alpha_Phi_minus_betaj_A_j`, and all
  `sqrtg_Ai` components.
- `ccc_BSSN` calls `ghl_interpolate_with_cell_centered_BSSN` through
  [`compute_ccc_BSSN.c`](../../../Unit_Tests/compute_ccc_BSSN.c). The unit test
  compares `alpha`, all shift components, `alpha_Phi_minus_betaj_A_j`, and all
  `sqrtg_Ai` components.
- `vvv_ADM` calls `ghl_interpolate_with_vertex_centered_ADM` through
  [`compute_vvv_ADM.c`](../../../Unit_Tests/compute_vvv_ADM.c). The unit test
  compares `alpha_Phi_minus_betaj_A_j` and all `sqrtg_Ai` components only.

No visible `vvv_BSSN` Induction interpolator test exists in
[`Unit_Tests/`](../../../Unit_Tests/) or
[`GRHayL/Induction/`](../../../GRHayL/Induction/).

## Repo-Local References

- [`GRHayL/include/ghl_induction.h`](../../../GRHayL/include/ghl_induction.h)
- [`GRHayL/Induction/Interpolators/interpolate_with_cell_centered_ADM.c`](../../../GRHayL/Induction/Interpolators/interpolate_with_cell_centered_ADM.c)
- [`GRHayL/Induction/Interpolators/interpolate_with_cell_centered_BSSN.c`](../../../GRHayL/Induction/Interpolators/interpolate_with_cell_centered_BSSN.c)
- [`GRHayL/Induction/Interpolators/interpolate_with_vertex_centered_ADM.c`](../../../GRHayL/Induction/Interpolators/interpolate_with_vertex_centered_ADM.c)
- [`GRHayL/Induction/Interpolators/interpolate_helper.c`](../../../GRHayL/Induction/Interpolators/interpolate_helper.c)
- [`GRHayL/Induction/Interpolators/ghl_induction_helpers.h`](../../../GRHayL/Induction/Interpolators/ghl_induction_helpers.h)
- [`GRHayL/Induction/Interpolators/make.code.defn`](../../../GRHayL/Induction/Interpolators/make.code.defn)
- [`Unit_Tests/compute_ccc_ADM.c`](../../../Unit_Tests/compute_ccc_ADM.c)
- [`Unit_Tests/compute_ccc_BSSN.c`](../../../Unit_Tests/compute_ccc_BSSN.c)
- [`Unit_Tests/compute_vvv_ADM.c`](../../../Unit_Tests/compute_vvv_ADM.c)
- [`Unit_Tests/unit_test_induction_ccc_ADM.c`](../../../Unit_Tests/unit_test_induction_ccc_ADM.c)
- [`Unit_Tests/unit_test_induction_ccc_BSSN.c`](../../../Unit_Tests/unit_test_induction_ccc_BSSN.c)
- [`Unit_Tests/unit_test_induction_vvv_ADM.c`](../../../Unit_Tests/unit_test_induction_vvv_ADM.c)
- [`docs/raw/Induction.dox`](../../../docs/raw/Induction.dox)
