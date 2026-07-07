# Induction Gem

## Purpose

Induction computes magnetic-vector-potential flux contributions, gauge terms, and densitized scalar-potential RHS terms for magnetic evolution.

## Read First

- `docs/raw/Induction.dox`
- `docs/raw/derivation.md`
- `GRHayL/include/ghl_induction.h`
- `wiki/physics/variables-and-conventions.md`
- `wiki/physics/evolution-equation-map.md`

## Public Headers

- `GRHayL/include/ghl_induction.h`
- `GRHayL/include/ghl.h`

Key public surface:
- `ghl_HLL_flux_with_B`
- `ghl_HLL_flux_with_Btilde`
- `ghl_interpolate_with_cell_centered_ADM`
- `ghl_interpolate_with_cell_centered_BSSN`
- `ghl_interpolate_with_vertex_centered_ADM`
- `ghl_calculate_phitilde_rhs`
- `ghl_HLL_vars`
- `ghl_induction_interp_vars`

## Implementation Paths

- Magnetic HLL flux: `GRHayL/Induction/HLL_flux_with_B.c`, `GRHayL/Induction/HLL_flux_with_Btilde.c`
- Scalar-potential RHS: `GRHayL/Induction/calculate_phitilde_rhs.c`
- Interpolators: `GRHayL/Induction/Interpolators/`
- Build lists: `GRHayL/Induction/make.code.defn`, `GRHayL/Induction/Interpolators/make.code.defn`

## Test Paths

- `Unit_Tests/unit_test_induction_ccc_ADM.c`
- `Unit_Tests/unit_test_induction_ccc_BSSN.c`
- `Unit_Tests/unit_test_induction_vvv_ADM.c`
- `Unit_Tests/unit_test_ET_Legacy_induction_gauge_rhs.c`
- `Unit_Tests/compute_A_flux_with_B.c`
- `Unit_Tests/compute_A_flux_with_Btilde.c`
- `Unit_Tests/compute_ccc_ADM.c`
- `Unit_Tests/compute_ccc_BSSN.c`
- `Unit_Tests/compute_vvv_ADM.c`
- `Unit_Tests/data_gen/unit_test_data_induction_interpolation.c`
- `Unit_Tests/data_gen/unit_test_data_ET_Legacy_induction_gauge_rhs.c`

## Key Contracts

- `A_i` is staggered in directions perpendicular to its component.
- `tilde{Phi}` is fully staggered.
- `B` may be staggered by direction; GRHayL only sees the stencils passed into its functions.
- `ghl_HLL_vars` uses cross-product component labels `1` and `2`; do not treat those names as fixed coordinate indices.
- `ghl_HLL_flux_with_B` needs `psi6`; `ghl_HLL_flux_with_Btilde` expects densitized magnetic input.
- `Lorenz_damping_factor` controls scalar-potential damping in `ghl_calculate_phitilde_rhs`.

## Common Edit Routes

- Change HLL flux formula: update both `B` and `Btilde` variants or explain why only one changes, then update flux tests.
- Change interpolation assumptions: update the matching interpolator, helper, Doxygen grid-location text, and interpolation tests.
- Change scalar-potential RHS: update `calculate_phitilde_rhs`, gauge RHS tests, and Induction Doxygen.
- Add metric staggering support: add interpolator, declarations, build list, tests, and docs.

## Drift Risks

- Staggering assumptions are easy to break because callers pass stencils without grid metadata.
- Characteristic-speed changes in Flux_Source can affect induction HLL inputs.
- Downstream GRHayLib includes `ghl_induction.h` and compiles Induction subdirectories; new source locations require coordination.

## Do Not Duplicate

Keep stencil diagrams, formulas, and function details in `docs/raw/Induction.dox` and source comments. This page stays routing-focused.
