# Induction HLL Flux

> Cross-product `ghl_HLL_vars` packing, `ghl_HLL_flux_with_B` versus `ghl_HLL_flux_with_Btilde`, stencil sample offsets, and unguarded speed denominators. · Status: confirmed
> Up: [Induction](index.md)

## Summary

`GRHayL/Induction/HLL_flux_with_B.c` and `HLL_flux_with_Btilde.c` compute scalar electric-field HLL flux contributions for staggered vector-potential components `A_x`, `A_y`, `A_z` in both default and `--disable-hdf5` builds. Both take caller-packed `ghl_HLL_vars` (using cyclic cross-product labels `1` and `2`), perform no input validation or denominator guards, and differ strictly by whether `psi6 = sqrt(gamma)` densitization is applied internally (`with_B`) or already present in magnetic inputs (`with_Btilde`).

## Detail

### `ghl_HLL_vars` Cyclic Component Packing And Densitization

| RHS component | Helper `A_dir` | Label `1` | Label `2` | Helper sample offsets (`v1*`/`v2*`, `B1*`, `B2*`) |
| --- | ---: | --- | --- | --- |
| `A_x` | `1` | `y` | `z` | `v`: `(i, j+1, k+1)`; `B1`: `(i, j, k+1)`; `B2`: `(i, j+1, k)` (Ev: `Unit_Tests/compute_A_flux_with_B.c`) |
| `A_y` | `2` | `z` | `x` | `v`: `(i+1, j, k+1)`; `B1`: `(i+1, j, k)`; `B2`: `(i, j, k+1)` (Ev: `Unit_Tests/compute_A_flux_with_Btilde.c`) |
| `A_z` | `3` | `x` | `y` | `v`: `(i+1, j+1, k)`; `B1`: `(i, j+1, k)`; `B2`: `(i+1, j, k)` |

In `compute_A_flux_with_{B,Btilde}.c` (which evaluate interior `2 <= i,j,k < dirlength - 2`), `c1_*` comes from label `1` at the `B2` offset and `c2_*` from label `2` at the `B1` offset. `ghl_HLL_flux_with_B` multiplies the complete scalar HLL result by caller-supplied `psi6 = sqrt(gamma)` at the `A_i` point (never densitizes individual packed face values), whereas `ghl_HLL_flux_with_Btilde` assumes magnetic fields are already densitized `Btilde = sqrt(gamma) * B`. Packed velocities are transport velocities `u^i / u^0`, and magnetic fields use GRHayL's `1 / sqrt(4*pi)` normalization without conversion.

Evidence (scientific/interface requirement): `GRHayL/include/ghl_induction.h::ghl_HLL_vars`; corroboration: `GRHayL/Induction/HLL_flux_with_B.c::ghl_HLL_flux_with_B` and `GRHayL/Induction/HLL_flux_with_Btilde.c::ghl_HLL_flux_with_Btilde`.

### Caller Hazards And Test Boundaries

- **Orientation & suffix order:** Swapping `B1`/`B2` flips cross-product sign. On velocity suffixes `rr`, `rl`, `lr`, `ll`, the first letter corresponds to label `1` and the second to label `2`; `B1` is reconstructed along label `2`, and `B2` along label `1`.
- **Unguarded denominators:** Both functions divide directly by `c1_min + c1_max` and `c2_min + c2_max` without checking for zero, flooring negative roundoff residues from `Flux_Source`, or guarding product overflow. Callers must ensure finite admissible bounds and usable denominators.
- **Test boundaries:** `unit_test_HLL_flux.c` exercises both `B` and `Btilde` across `Ax_rhs`, `Ay_rhs`, `Az_rhs`, while `unit_test_ET_Legacy_HLL_flux.c` covers only `B` against external legacy data. Because `unit_test_data_HLL_flux.c` draws random signed speeds in `[-10, 10]` (rejecting only exact zero sums), it tests HLL algebra/packing rather than Flux_Source speed coupling or zero-denominator safety. No checked-in call site appears in `implementations/GRHayLib/`.

## Sources

- [GRHayL/include/ghl_induction.h](../../../GRHayL/include/ghl_induction.h) - `ghl_HLL_vars`, `ghl_HLL_flux_with_B`, `ghl_HLL_flux_with_Btilde`
- [GRHayL/Induction/HLL_flux_with_B.c](../../../GRHayL/Induction/HLL_flux_with_B.c) - undensitized `B` + `psi6` flux helper
- [GRHayL/Induction/HLL_flux_with_Btilde.c](../../../GRHayL/Induction/HLL_flux_with_Btilde.c) - densitized `Btilde` flux helper
- [GRHayL/Induction/make.code.defn](../../../GRHayL/Induction/make.code.defn) - Induction build manifest
- [Unit_Tests/compute_A_flux_with_B.c](../../../Unit_Tests/compute_A_flux_with_B.c) - `A_dir`, `v_offset`, `B1_offset`, `B2_offset`, `psi6` packing
- [Unit_Tests/compute_A_flux_with_Btilde.c](../../../Unit_Tests/compute_A_flux_with_Btilde.c) - `Btilde` packing helper
- [Unit_Tests/unit_test_HLL_flux.c](../../../Unit_Tests/unit_test_HLL_flux.c) - local `B` and `Btilde` HLL test
- [Unit_Tests/unit_test_ET_Legacy_HLL_flux.c](../../../Unit_Tests/unit_test_ET_Legacy_HLL_flux.c) - ET_Legacy `B` HLL test
- [docs/raw/Induction.dox](../../../docs/raw/Induction.dox) - Induction equations and staggering docs

## See Also

- Depends on: [Characteristic Speeds](../flux-source/characteristic-speeds.md)
- Depends on: [Variables And Conventions](../../physics/variables-and-conventions.md)
- Validated by: [Induction Tests And Reference Data](tests-and-reference-data.md)
- See also: [Interpolation And Staggering](interpolation-and-staggering.md)
