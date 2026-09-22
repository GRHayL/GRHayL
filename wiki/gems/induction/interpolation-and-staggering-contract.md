# Induction Interpolation And Staggering Contract

## Routing Purpose

Use this page when checking Induction interpolator API shape,
vector-potential staggering, caller packing, output fields, downstream usage,
or test coverage. Start from the [Induction hub](../induction.md), then verify
the public signatures in
[`ghl_induction.h`](../../../GRHayL/include/ghl_induction.h) and behavior in
the [interpolator source](../../../GRHayL/Induction/Interpolators/).

GRHayL receives plain caller-packed arrays. It has no grid metadata from which
to infer centering, index origin, or ghost-zone reach. The caller must select
the matching wrapper and pack its exact stencil.

## Public API

The public interpolators are:

- `ghl_interpolate_with_cell_centered_ADM`: [source](../../../GRHayL/Induction/Interpolators/interpolate_with_cell_centered_ADM.c),
  [test helper](../../../Unit_Tests/compute_ccc_ADM.c), and
  [test](../../../Unit_Tests/unit_test_induction_ccc_ADM.c).
- `ghl_interpolate_with_cell_centered_BSSN`: [source](../../../GRHayL/Induction/Interpolators/interpolate_with_cell_centered_BSSN.c),
  [test helper](../../../Unit_Tests/compute_ccc_BSSN.c), and
  [test](../../../Unit_Tests/unit_test_induction_ccc_BSSN.c).
- `ghl_interpolate_with_vertex_centered_ADM`: backward/current
  vertex-centered ADM interpolation;
  [source](../../../GRHayL/Induction/Interpolators/interpolate_with_vertex_centered_ADM.c),
  [fixture helper](../../../Unit_Tests/compute_vvv_ADM.c), and
  [test](../../../Unit_Tests/unit_test_induction_vvv_ADM.c).

The output type `ghl_induction_interp_vars` contains `alpha`,
`alpha_Phi_minus_betaj_A_j`, `sqrtg_Ai[3]`, and `betai[3]`. Cell-centered
ADM and BSSN wrappers assign every field. The vertex-centered wrapper assigns
only `sqrtg_Ai[3]` and `alpha_Phi_minus_betaj_A_j`; `alpha` and `betai[3]`
retain their incoming or uninitialized values. Callers must obtain lapse and
shift directly from the current vertex when computing a later `phitilde` RHS.

None of these routines returns an error, validates bounds, or checks centering.

## Shared Array Layout

C array order is `[z][y][x]`. The public signatures use:

| Input | Shape | Use |
| --- | --- | --- |
| `metric_stencil` or `metric` | `2x2x2` `ghl_metric_quantities` | every wrapper |
| `psi_stencil` | `2x2x2` scalar | cell-centered BSSN only |
| `Ax_stencil` | `3x3x3` scalar | every wrapper |
| `Ay_stencil` | `3x3x3` scalar | every wrapper |
| `Az_stencil` | `3x3x3` scalar | every wrapper |
| `phitilde` | scalar | every wrapper |

Relative to `tildePhi(i+1/2,j+1/2,k+1/2)`, vector-potential arrays cover:

| Array | x range | y range | z range |
| --- | --- | --- | --- |
| `Ax_stencil[z][y][x]` | `i-1..i+1` | `j-1/2..j+3/2` | `k-1/2..k+3/2` |
| `Ay_stencil[z][y][x]` | `i-1/2..i+3/2` | `j-1..j+1` | `k-1/2..k+3/2` |
| `Az_stencil[z][y][x]` | `i-1/2..i+3/2` | `j-1/2..j+3/2` | `k-1..k+1` |

ADM stencil elements must have `lapse`, `betaU`, `gammaUU`, and
`sqrt_detgamma` initialized. BSSN callers instead provide conformal
inverse-metric components in `metric_stencil.gammaUU` and physical `psi`
separately.

## Vertex-Centered ADM Contract

`ghl_interpolate_with_vertex_centered_ADM` uses this metric packing:

| Role | Array entry | Physical vertex |
| --- | --- | --- |
| current `tildePhi` vertex | `[1][1][1]` | `(i+1/2,j+1/2,k+1/2)` |
| backward x neighbor | `[1][1][0]` | `(i-1/2,j+1/2,k+1/2)` |
| backward y neighbor | `[1][0][1]` | `(i+1/2,j-1/2,k+1/2)` |
| backward z neighbor | `[0][1][1]` | `(i+1/2,j+1/2,k-1/2)` |

The caller therefore packs offsets `-1..0` relative to the current
`tildePhi` vertex on each axis. In the unit-test helper's array index space,
these are indices `[i-1..i]`, `[j-1..j]`, and `[k-1..k]`. Only the four
listed entries are read. The fixed `2x2x2` parameter shape does not mean all
eight entries contribute.

`A_x`, `A_y`, and `A_z` each live one half-step behind `tildePhi` along their
component axis. Raising the index at an `A_i` location therefore needs the
metric average between the current `tildePhi` vertex and the backward vertex
on that axis. The three pairs are:

- x: `[1][1][0]` and `[1][1][1]`;
- y: `[1][0][1]` and `[1][1][1]`;
- z: `[0][1][1]` and `[1][1][1]`.

Using the current vertex and a forward neighbor instead moves each metric
average one complete grid cell ahead of its corresponding `A_i` location.
That is a numerical displacement, not an array-order convention. The public
function therefore has one contract: backward/current packing with the current
vertex at `[1][1][1]`.

## Downstream Evidence And Limits

The original
[WVUThorns IllinoisGRMHD gauge implementation](https://bitbucket.org/zach_etienne/wvuthorns/src/master/IllinoisGRMHD/src/Lorenz_psi6phi_rhs__add_gauge_terms_to_A_i_rhs.C)
implements its interpolation directly. It does not call the GRHayL vertex
wrapper, so it is historical evidence for the algorithm family rather than
authority for this function's argument contract.

[GRHayLET `main` IllinoisGRMHD](https://github.com/GRHayL/GRHayLET/blob/main/IllinoisGRMHD/src/evaluate_phitilde_and_A_gauge_rhs.c)
packs a cell-centered `2x2x2` metric stencil and calls
`ghl_interpolate_with_cell_centered_ADM`. It does not consume the vertex API.
The [GRHayLHD](https://github.com/GRHayL/GRHayLET/tree/main/GRHayLHD) and
[GRHayLHDX](https://github.com/GRHayL/GRHayLET/tree/main/GRHayLHDX) thorns are
hydrodynamic integrations; no call to the vertex interpolator is visible in
those trees.

The experimental
[IllinoisGRMHDX branch](https://github.com/GRHayL/GRHayLET/blob/scupp/GRHayLMHDX/IllinoisGRMHDX/src/evaluate_phitilde_and_A_gauge_rhs.cxx)
packs the metric over `-1..0`, puts the current vertex at `[1][1][1]`, and
calls `ghl_interpolate_with_vertex_centered_ADM`. Its packing matches this
contract. The same file leaves gauge derivatives disabled and sets
`phitilde_rhs` to zero, so this development branch is evidence of intended
centering, not completed runtime validation.

Searches of public GRHayL and GRHayLET trees, GRHayLHD, GRHayLHDX, pyghl,
and visible GRHayL forks found no other direct downstream call. This bounds
public evidence; it cannot rule out private or unindexed consumers.

## Internal Helpers

Wrappers call private routines declared in
[`ghl_induction_helpers.h`](../../../GRHayL/Induction/Interpolators/ghl_induction_helpers.h)
and implemented in
[`interpolate_helper.c`](../../../GRHayL/Induction/Interpolators/interpolate_helper.c):

- `ghl_A_i_avg` averages each vector-potential component to `tildePhi`,
  `A_x`, `A_y`, and `A_z` locations.
- `ghl_ADM_cell_interp` handles cell-centered ADM interpolation.
- `ghl_BSSN_cell_interp` handles cell-centered BSSN and `psi` interpolation.
- `ghl_ADM_vertex_interp` averages the three backward/current axial metric
  pairs listed above.

These are internal machinery. External callers should use public wrappers.

## Coverage

- The cell-centered ADM and BSSN fixture replays compare lapse, shift,
  `alpha_Phi_minus_betaj_A_j`, and every `sqrtg_Ai` component.
- The `vvv_ADM` fixture replay packs metric offsets `-1..0` and calls
  `ghl_interpolate_with_vertex_centered_ADM`. The fixture is a regression
  reference generated through the same helper and wrapper path, not an
  independent numerical oracle.
- The direct asymmetric affine case calls the same public function. It checks
  every `sqrtg_Ai` component, the gauge scalar, coordinate orientation, and
  that `alpha` and `betai` remain unchanged.

Test helpers fill complete `3x3x3` vector-potential arrays and run only for
`1 <= i,j,k < dirlength-1`. Cell-centered metric helpers pack offsets
`0..+1`; the vertex-centered helper packs `-1..0`.

No public `vvv_BSSN` declaration, implementation, helper, or test is visible.

## Repo-Local References

- [Public Induction header](../../../GRHayL/include/ghl_induction.h)
- [Vertex ADM wrapper](../../../GRHayL/Induction/Interpolators/interpolate_with_vertex_centered_ADM.c)
- [Interpolation helpers](../../../GRHayL/Induction/Interpolators/interpolate_helper.c)
- [Internal helper declarations](../../../GRHayL/Induction/Interpolators/ghl_induction_helpers.h)
- [Vertex fixture helper](../../../Unit_Tests/compute_vvv_ADM.c)
- [Vertex test](../../../Unit_Tests/unit_test_induction_vvv_ADM.c)
- [Induction Doxygen source](../../../docs/raw/Induction.dox)
- [Induction tests and fixtures](tests-and-fixtures.md)
