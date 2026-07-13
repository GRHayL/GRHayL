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
from that wrapper. Those fields retain prior/uninitialized storage. None of the
wrappers returns an error or checks stencil bounds. The wrappers do not infer
grid locations from metadata.

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

C array order is `[z][y][x]`, matching test packing. Relative to the fully
staggered `tildePhi` point `(i+1/2,j+1/2,k+1/2)`, physical ranges are:

| Array | x range | y range | z range |
| --- | --- | --- | --- |
| `Ax_stencil[z][y][x]` | `i-1..i+1` | `j-1/2..j+3/2` | `k-1/2..k+3/2` |
| `Ay_stencil[z][y][x]` | `i-1/2..i+3/2` | `j-1..j+1` | `k-1/2..k+3/2` |
| `Az_stencil[z][y][x]` | `i-1/2..i+3/2` | `j-1/2..j+3/2` | `k-1..k+1` |

Cell-centered metric and BSSN `psi` stencils span `[i..i+1]`, `[j..j+1]`,
and `[k..k+1]`; all eight entries contribute to lapse/shift quantities. The
vertex-centered ADM signature is also `2x2x2`, but implementation reads only
`[0][0][0]`, `[0][0][1]`, `[0][1][0]`, and `[1][0][0]`: current `tildePhi`
vertex plus forward x/y/z neighbors.

`A_i` is staggered in directions perpendicular to its component, and
`tildePhi` is fully staggered. An opening comment in `interpolate_helper.c`
mistypes the `A_z` location as `(i+1/2,j,k+1/2)`; function Doxygen, wrapper
comments, array accesses, and tests consistently use `(i+1/2,j+1/2,k)`.

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

ADM stencil elements must have `lapse`, `betaU`, `gammaUU`, and
`sqrt_detgamma` initialized as used by each wrapper. BSSN callers provide
conformal inverse-metric components in `metric_stencil.gammaUU` plus a separate
physical `psi_stencil`; this wrapper is not accepting a complete physical ADM
metric despite reusing `ghl_metric_quantities` as its container.

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

More precisely, no public `vvv_BSSN` declaration, definition, manifest entry,
helper, or test exists. This is an unsupported variant, not merely an untested
existing API.

Test helpers fill complete `3x3x3` A stencils and complete cell-centered metric
stencils, then call wrappers only for `1 <= i,j,k < dirlength-1`; this bounds
their `-1..+1` A reach and `0..+1` metric reach. Fixture outputs are produced by
the same helper/wrapper path and are regression references rather than an
independent interpolation oracle. Runner and all compiler workflows configure
the three supported replays; tracked files alone do not establish a current
execution result.

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
