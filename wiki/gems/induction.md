# Induction Gem

## Purpose

Induction computes magnetic-vector-potential flux contributions, gauge terms, and densitized scalar-potential RHS terms for magnetic evolution.

## Read First

- `docs/raw/Induction.dox` as read-only evidence for KB work
- `docs/raw/derivation.md`
- `GRHayL/include/ghl_induction.h`
- `wiki/physics/variables-and-conventions.md`
- `wiki/physics/evolution-equation-map.md`

## KB Routes

- [HLL flux contract](induction/hll-flux-contract.md)
- [Interpolation and staggering contract](induction/interpolation-and-staggering-contract.md)
- [Gauge RHS contract](induction/gauge-rhs-contract.md)
- [Tests and fixtures](induction/tests-and-fixtures.md)
- [Verification workflows](induction/verification-workflows.md)

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

All direct Induction routines are declared, defined, and in normal/no-HDF5
manifests. GRHayLib lists their directories for compilation but contains no
checked-in call sites. Tests and workflow routes are configured evidence, not
an observed pass.

## Implementation Paths

- Magnetic HLL flux: `GRHayL/Induction/HLL_flux_with_B.c`, `GRHayL/Induction/HLL_flux_with_Btilde.c`
- Scalar-potential RHS: `GRHayL/Induction/calculate_phitilde_rhs.c`
- Interpolators: `GRHayL/Induction/Interpolators/`
- Build lists: `GRHayL/Induction/make.code.defn`, `GRHayL/Induction/Interpolators/make.code.defn`

## Test Paths

- `Unit_Tests/unit_test_HLL_flux.c`
- `Unit_Tests/unit_test_ET_Legacy_HLL_flux.c`
- `Unit_Tests/unit_test_induction_ccc_ADM.c`
- `Unit_Tests/unit_test_induction_ccc_BSSN.c`
- `Unit_Tests/unit_test_induction_vvv_ADM.c`
- `Unit_Tests/unit_test_ET_Legacy_induction_gauge_rhs.c`
- `Unit_Tests/compute_A_flux_with_B.c`
- `Unit_Tests/compute_A_flux_with_Btilde.c`
- `Unit_Tests/compute_ccc_ADM.c`
- `Unit_Tests/compute_ccc_BSSN.c`
- `Unit_Tests/compute_vvv_ADM.c`
- `Unit_Tests/data_gen/unit_test_data_HLL_flux.c`
- `Unit_Tests/data_gen/unit_test_data_ET_Legacy_HLL_flux.c`
- `Unit_Tests/data_gen/unit_test_data_induction_interpolation.c`
- `Unit_Tests/data_gen/unit_test_data_ET_Legacy_induction_gauge_rhs.c`

## Key Contracts

- `A_i` is staggered in directions perpendicular to its component.
- `tilde{Phi}` is fully staggered.
- `B` may be staggered by direction; GRHayL only sees the stencils passed into its functions.
- `ghl_HLL_vars` uses cross-product component labels `1` and `2`; do not treat those names as fixed coordinate indices.
- `ghl_HLL_flux_with_B` needs `psi6`; `ghl_HLL_flux_with_Btilde` expects densitized magnetic input.
- HLL denominators require nonzero `cmin+cmax`; visible HLL fixtures use random
  speeds and do not test Flux_Source-to-Induction coupling.
- `Lorenz_damping_factor` controls scalar-potential damping in `ghl_calculate_phitilde_rhs`.
- Vertex-centered ADM interpolation intentionally leaves output `alpha` and
  `betai` unwritten; see interpolation contract before unpacking output.

## Common Edit Routes

- Change HLL flux formula: use [HLL flux contract](induction/hll-flux-contract.md), check both `B` and `Btilde` variants or explain why only one changes, then route tests through [tests and fixtures](induction/tests-and-fixtures.md).
- Change interpolation assumptions: use [interpolation and staggering contract](induction/interpolation-and-staggering-contract.md), then check matching interpolator, helper, declaration, build list, and interpolation tests.
- Change scalar-potential RHS: use [gauge RHS contract](induction/gauge-rhs-contract.md), then check `calculate_phitilde_rhs`, ET Legacy gauge RHS coverage, and targeted runs in [verification workflows](induction/verification-workflows.md).
- Add metric staggering support: check public declarations, interpolator source, local build lists, and test/fixture routing.

## Drift Risks

- Staggering assumptions are easy to break because callers pass stencils without grid metadata.
- Characteristic-speed changes in Flux_Source can affect induction HLL inputs.
- Downstream GRHayLib includes `ghl_induction.h` and compiles Induction subdirectories; new source locations require coordination.
- Legacy `generate_makefile.sh` emits bogus Induction targets; normal
  `configure` parsing remains separate. Use [verification workflows](induction/verification-workflows.md).

## Do Not Duplicate

Keep stencil diagrams, formulas, and function details in source/Doxygen evidence. This page stays routing-focused.
