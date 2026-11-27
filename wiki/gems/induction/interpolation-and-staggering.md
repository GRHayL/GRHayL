# Induction Interpolation And Staggering

> Stencil packing (`2x2x2` metric, `3x3x3` `A_i`), output fields (`ghl_induction_interp_vars`), backward/current vertex-centered ADM indexing, and downstream GRHayLET alignment. · Status: confirmed
> Up: [Induction](index.md)

## Summary

`GRHayL/Induction/Interpolators/` provides public staggered metric/vector-potential interpolators (`ghl_interpolate_with_cell_centered_ADM`, `ghl_interpolate_with_cell_centered_BSSN`, `ghl_interpolate_with_vertex_centered_ADM`) receiving caller-packed `[z][y][x]` arrays without grid metadata or validation. Cell-centered wrappers write every field of `ghl_induction_interp_vars`; vertex-centered ADM writes only `sqrtg_Ai[3]` and `alpha_Phi_minus_betaj_A_j`, leaving `alpha` and `betai[3]` untouched.

## Detail

### Public Wrappers, Output Fields, And Array Packing

`ghl_induction_interp_vars` (`GRHayL/include/ghl_induction.h`) holds `alpha`, `betai[3]`, `alpha_Phi_minus_betaj_A_j` (at `phitilde = sqrt(gamma) * Phi`), and `sqrtg_Ai[3]` (representing `alpha * sqrt(gamma) * A^i = sqrt(-g) * A^i` at each staggered `A_i`, lapse already included). Each public wrapper returns `void` without validation (ADM divides by `sqrt_detgamma`; BSSN divides by `psi^6`). No `vvv_BSSN` wrapper exists.

- **Array order (`[z][y][x]`):** `metric_stencil` (`2x2x2` `ghl_metric_quantities`), `psi_stencil` (`2x2x2` for cell-centered BSSN only), scalar `phitilde`, and `3x3x3` `Ax_stencil` (`i-1..i+1`, `j-1/2..j+3/2`, `k-1/2..k+3/2`), `Ay_stencil` (`i-1/2..i+3/2`, `j-1..j+1`, `k-1/2..k+3/2`), `Az_stencil` (`i-1/2..i+3/2`, `j-1/2..j+3/2`, `k-1..k+1`) relative to `tildePhi(i+1/2, j+1/2, k+1/2)`.
- **ADM vs. BSSN discretization:** ADM stencils initialize `lapse`, `betaU`, `gammaUU`, `sqrt_detgamma` and average the combined `alpha * sqrt(gamma) * gamma^{ij}`. BSSN expects conformal inverse `psi^4 * gamma^{ij}` in `metric_stencil.gammaUU` (upper triangle only; ignores `sqrt_detgamma`), assumes a unit-determinant conformal metric (`psi^6` without `sqrt(det(gamma_tilde))`), and averages `alpha * psi^2` and conformal `gamma^{ij}` separately before multiplying—so ADM and BSSN need not match numerically on spatially varying geometry.
- **Internal helpers (`ghl_induction_helpers.h`, `interpolate_helper.c`):** `ghl_A_i_avg`, `ghl_ADM_cell_interp`, `ghl_BSSN_cell_interp`, and `ghl_ADM_vertex_interp` are internal machinery behind the public wrappers.

Evidence (descriptive behavior): `GRHayL/Induction/Interpolators/interpolate_with_vertex_centered_ADM.c::ghl_interpolate_with_vertex_centered_ADM`; corroboration: `GRHayL/Induction/Interpolators/interpolate_helper.c::ghl_ADM_vertex_interp`.

### Backward/Current Vertex-Centered ADM Packing And Downstream Evidence

Because `A_x`, `A_y`, `A_z` each sit half a cell behind `tildePhi` along their component axis, `ghl_interpolate_with_vertex_centered_ADM` reads **only four entries** of `metric_stencil[2][2][2]`, packing offsets `-1..0` with the current `tildePhi` vertex at `[1][1][1]`:
- `[1][1][1]` = `(i+1/2, j+1/2, k+1/2)` (current `tildePhi` vertex);
- `[1][1][0]` = backward-x `(i-1/2, j+1/2, k+1/2)` (averaged with `[1][1][1]` for `A_x`);
- `[1][0][1]` = backward-y `(i+1/2, j-1/2, k+1/2)` (averaged with `[1][1][1]` for `A_y`);
- `[0][1][1]` = backward-z `(i+1/2, j+1/2, k-1/2)` (averaged with `[1][1][1]` for `A_z`).

Packing forward offsets `0..+1` displaces every metric average by one full grid cell. Callers must obtain vertex `alpha` and `betai` directly when forming `phitilde` RHS. `unit_test_induction_vvv_ADM.c` includes a direct asymmetric affine check confirming all three `sqrtg_Ai` components, the gauge scalar, orientation, and unchanged `alpha`/`betai`.

**Downstream context:** `GRHayLET` `main/IllinoisGRMHD` calls `ghl_interpolate_with_cell_centered_ADM`; `GRHayLHD` and `GRHayLHDX` are hydro-only; the experimental `IllinoisGRMHDX` branch packs `-1..0` with `[1][1][1]` into `ghl_interpolate_with_vertex_centered_ADM` (while leaving gauge derivatives disabled and `phitilde_rhs = 0`); legacy `WVUThorns/IllinoisGRMHD` interpolates inline.

## Sources

- [GRHayL/include/ghl_induction.h](../../../GRHayL/include/ghl_induction.h) - public interpolators and `ghl_induction_interp_vars`
- [GRHayL/Induction/Interpolators/interpolate_with_cell_centered_ADM.c](../../../GRHayL/Induction/Interpolators/interpolate_with_cell_centered_ADM.c) - cell-centered ADM wrapper
- [GRHayL/Induction/Interpolators/interpolate_with_cell_centered_BSSN.c](../../../GRHayL/Induction/Interpolators/interpolate_with_cell_centered_BSSN.c) - cell-centered BSSN wrapper
- [GRHayL/Induction/Interpolators/interpolate_with_vertex_centered_ADM.c](../../../GRHayL/Induction/Interpolators/interpolate_with_vertex_centered_ADM.c) - backward/current vertex-centered ADM wrapper
- [GRHayL/Induction/Interpolators/interpolate_helper.c](../../../GRHayL/Induction/Interpolators/interpolate_helper.c) - `ghl_A_i_avg`, `ghl_ADM_cell_interp`, `ghl_BSSN_cell_interp`, `ghl_ADM_vertex_interp`
- [Unit_Tests/unit_test_induction_vvv_ADM.c](../../../Unit_Tests/unit_test_induction_vvv_ADM.c) - vertex-centered ADM replay and asymmetric affine check
- [WVUThorns IllinoisGRMHD gauge RHS](https://bitbucket.org/zach_etienne/wvuthorns/src/master/IllinoisGRMHD/src/Lorenz_psi6phi_rhs__add_gauge_terms_to_A_i_rhs.C) - historical algorithm family
- [GRHayLET main IllinoisGRMHD gauge RHS](https://github.com/GRHayL/GRHayLET/blob/main/IllinoisGRMHD/src/evaluate_phitilde_and_A_gauge_rhs.c) - cell-centered ADM caller
- [GRHayLET IllinoisGRMHDX gauge RHS](https://github.com/GRHayL/GRHayLET/blob/scupp/GRHayLMHDX/IllinoisGRMHDX/src/evaluate_phitilde_and_A_gauge_rhs.cxx) - experimental vertex-centered caller
- [GRHayLET GRHayLHD](https://github.com/GRHayL/GRHayLET/tree/main/GRHayLHD) - hydro-only thorn context
- [GRHayLET GRHayLHDX](https://github.com/GRHayL/GRHayLET/tree/main/GRHayLHDX) - CarpetX hydro-only thorn context

## See Also

- Depends on: [Gauge RHS](gauge-rhs.md)
- Validated by: [Induction Tests And Reference Data](tests-and-reference-data.md)
- See also: [HLL Flux](hll-flux.md)
