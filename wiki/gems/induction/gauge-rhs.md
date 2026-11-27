# Induction Gauge RHS

> Caller stencil layout, shift-upwinding branches, `Lorenz_damping_factor` handling, and ET_Legacy coverage for `ghl_calculate_phitilde_rhs`. · Status: confirmed
> Up: [Induction](index.md)

## Summary

`ghl_calculate_phitilde_rhs` (`GRHayL/Induction/calculate_phitilde_rhs.c`) returns the scalar-potential evolution RHS from caller-assembled 5-point shift/`phitilde` lines, 2-point forward `sqrtg_Ai` stencils, inverse spacings `dxi[3]`, and explicit `Lorenz_damping_factor`. There is no separate public function for the `A_i` gauge RHS contribution (callers difference `alpha_Phi_minus_betaj_A_j` directly) and no standalone `phitilde_rhs`-only unit test.

## Detail

### Arguments, Shift Upwinding, And Stencil Reach

Declared in `GRHayL/include/ghl_induction.h` and built in both default and `--disable-hdf5` modes (`GRHayL/Induction/make.code.defn`), `ghl_calculate_phitilde_rhs` performs no pointer or ghost-zone validation:
- `dxi[3]`: inverse grid spacings `1/dx`, `1/dy`, `1/dz` (not `dx`).
- `Lorenz_damping_factor`: passed directly as a scalar (not read from `ghl_parameters`; callers keep it consistent with `params->Lorenz_damping_factor`).
- `alpha_interp`: lapse at the `phitilde` location.
- `shift{x,y,z}_interp[5]` and `phitilde_stencil[3][5]`: 5-point directional lines spanning `-2..+2` around center index `[2]`. Upwinding branches on center shift `[2]`: negative center shift selects the backward-biased stencil, non-negative selects the forward-biased stencil. The damping term multiplies `Lorenz_damping_factor`, `alpha_interp`, and `phitilde_stencil[0][2]` (assuming all three directional lines share the same center point).
- `sqrtg_Ai_stencil[3][2]`: component-major `[component][0]` (current point) and `[component][1]` (forward neighbor along that component axis) for `sqrtg_Ai = alpha * sqrt(gamma) * A^i` (lapse already included).

Routine-local stencil reach is two points backward/forward along shift and `phitilde` lines and one forward point for `sqrtg_Ai`. Cell-centered interpolators supply all needed fields, whereas vertex-centered ADM callers supply vertex `alpha` and shift separately ([Interpolation And Staggering](interpolation-and-staggering.md)).

Evidence (descriptive behavior): `GRHayL/Induction/calculate_phitilde_rhs.c::ghl_calculate_phitilde_rhs`; corroboration: `GRHayL/include/ghl_induction.h::ghl_calculate_phitilde_rhs`.

### ET_Legacy Coverage And Gaps

`Unit_Tests/unit_test_ET_Legacy_induction_gauge_rhs.c` (input generator `unit_test_data_ET_Legacy_induction_gauge_rhs.c`) runs cell-centered BSSN interpolation on `1 <= i,j,k < dirlength - 1`, then assembles `A_i` gauge finite differences and calls `ghl_calculate_phitilde_rhs` with fixed damping `0.1` on `3 <= i,j,k < dirlength - 3` against externally supplied ET_Legacy outputs. No standalone `phitilde_rhs`-only unit test, ADM/vertex-centered gauge RHS test, or explicit center-mismatch test exists.

## Sources

- [GRHayL/include/ghl_induction.h](../../../GRHayL/include/ghl_induction.h) - `ghl_calculate_phitilde_rhs` declaration
- [GRHayL/Induction/calculate_phitilde_rhs.c](../../../GRHayL/Induction/calculate_phitilde_rhs.c) - shift upwinding, divergence, and damping formula
- [Unit_Tests/unit_test_ET_Legacy_induction_gauge_rhs.c](../../../Unit_Tests/unit_test_ET_Legacy_induction_gauge_rhs.c) - combined BSSN + gauge RHS replay
- [Unit_Tests/data_gen/unit_test_data_ET_Legacy_induction_gauge_rhs.c](../../../Unit_Tests/data_gen/unit_test_data_ET_Legacy_induction_gauge_rhs.c) - ET_Legacy input generator
- [docs/raw/Induction.dox](../../../docs/raw/Induction.dox) - Lorenz gauge equations and caller-side `A_i` term

## See Also

- Depends on: [Interpolation And Staggering](interpolation-and-staggering.md)
- Validated by: [Induction Tests And Reference Data](tests-and-reference-data.md)
- See also: [Induction Verification Workflows](verification-workflows.md)
