# PPM Flow

> Piecewise parabolic reconstruction wrappers, 5-point helpers, shock flattening (`ftilde`), contact steepening, `ghl_parameters` controls, and ET_Legacy test boundaries. · Status: confirmed
> Up: [Reconstruction](index.md)

## Summary

`GRHayL/Reconstruction/PPM/` implements 6-point left-face PPM wrappers (`ghl_ppm_reconstruction`, `ghl_ppm_reconstruction_with_steepening`) built from two shifted 5-point centered helper calls, paired with `ghl_compute_ftilde` shock flattening and optional density contact steepening controlled by seven `ppm_*` fields in `ghl_parameters`. No dedicated standalone PPM unit test exists; current test coverage is through `unit_test_ET_Legacy_reconstruction.c`.

## Detail

### Wrapper And Helper Call Sequence

- **`ghl_ppm_reconstruction` (`ppm_reconstruction.c`):** accepts `ftilde[2]` and `var_data[6]` (`i-3 .. i+2`), calls `ghl_ppm_compute_for_cell(ftilde[0], var_data, ...)` taking helper right output into `var_datal`, then calls `ghl_ppm_compute_for_cell(ftilde[1], &var_data[1], ...)` taking shifted helper left output into `var_datar`.
- **`ghl_ppm_reconstruction_with_steepening` (`ppm_reconstruction_with_steepening.c`):** follows the same two-call structure while passing `params`, aligned 5-point `pressure` slices, and caller-supplied `Gamma_eff` to `ghl_ppm_compute_for_cell_with_steepening`.
- **5-point helpers (`ppm_compute_for_cell.c`, `ppm_compute_for_cell_with_steepening.c`, `slope_limit.c`):** compute limited centered slopes via `ghl_slope_limit` (zero when adjacent differences change sign), form provisional parabolic face values, optionally call `ghl_steepen_var`, apply flattening with caller-supplied `ftilde`, and enforce monotonicity.

Evidence (descriptive behavior): `GRHayL/Reconstruction/PPM/ppm_reconstruction.c::ghl_ppm_reconstruction`; corroboration: `GRHayL/Reconstruction/PPM/ppm_compute_for_cell_with_steepening.c::ghl_ppm_compute_for_cell_with_steepening`.

### Flattening, Steepening, Parameters, And Coverage Limits

- **Flattening (`ppm_compute_ftilde.c`, `shock_detection_ftilde.c`):** `ghl_compute_ftilde` evaluates `ghl_shock_detection_ftilde` on adjacent 5-point slices of `pressure[6]` and `v_flux_dirn[6]` to fill `ftilde[2]`. Shock detection compares `fabs(dP) / average` against `1e-15` for roundoff suppression and reads `ppm_flattening_epsilon`, `ppm_flattening_omega1`, and `ppm_flattening_omega2`. Centered PPM helpers apply `ftilde` literally without clamping—callers must supply finite values in `[0, 1]` (`0` leaves provisional values unflattened; `1` collapses to cell center before monotonization).
- **Contact steepening (`steepen_var.c`):** intended for density (pressure and velocity use the non-steepened wrapper in ET_Legacy). `ghl_steepen_var` reads `ppm_shock_k0`, `ppm_shock_eta1`, `ppm_shock_eta2`, `ppm_shock_epsilon`, and caller-provided `Gamma_eff` (no EOS lookup). When pressure/density minima are positive and ratios are finite, it uses normalized jump ratios; otherwise it falls back to a product-form inequality without rejecting invalid inputs.
- **Parameters & downstream coupling:** Core `initialize_params.c` sets defaults (`0.33`, `0.75`, `10.0`; `k0 = 0.1`, `eta1 = 20.0`, `eta2 = 0.05`, `epsilon = 0.01`), which `implementations/GRHayLib/` exposes and overwrites from `param.ccl` ([Struct Pack And Unpack](../../core/struct-pack-unpack.md), [GRHayLib Runtime Parameters](../../implementations/grhaylib-runtime-parameters.md)). PPM routines do not validate `ghl_parameters`.
- **Coverage boundary:** `unit_test_ET_Legacy_reconstruction.c` (using default PPM parameters and polytropic/thermal `Gamma = 2`) exercises both wrappers and `ghl_compute_ftilde` transitively; there is no dedicated `unit_test_PPM_reconstruction.c`, no PPM data generator, and no isolated test of the 5-point helpers or varying `Gamma_eff`.

Evidence (descriptive behavior): `GRHayL/Reconstruction/PPM/shock_detection_ftilde.c::ghl_shock_detection_ftilde`; corroboration: `GRHayL/Reconstruction/PPM/steepen_var.c::ghl_steepen_var`.

## Sources

- [GRHayL/include/ghl_reconstruction.h](../../../GRHayL/include/ghl_reconstruction.h) - PPM declarations
- [GRHayL/Reconstruction/PPM/make.code.defn](../../../GRHayL/Reconstruction/PPM/make.code.defn) - PPM build manifest
- [GRHayL/Reconstruction/PPM/ppm_reconstruction.c](../../../GRHayL/Reconstruction/PPM/ppm_reconstruction.c) - non-steepened 6-point wrapper
- [GRHayL/Reconstruction/PPM/ppm_reconstruction_with_steepening.c](../../../GRHayL/Reconstruction/PPM/ppm_reconstruction_with_steepening.c) - steepened 6-point wrapper
- [GRHayL/Reconstruction/PPM/ppm_compute_for_cell.c](../../../GRHayL/Reconstruction/PPM/ppm_compute_for_cell.c) - centered 5-point helper
- [GRHayL/Reconstruction/PPM/ppm_compute_for_cell_with_steepening.c](../../../GRHayL/Reconstruction/PPM/ppm_compute_for_cell_with_steepening.c) - steepened 5-point helper
- [GRHayL/Reconstruction/PPM/ppm_compute_ftilde.c](../../../GRHayL/Reconstruction/PPM/ppm_compute_ftilde.c) - `ghl_compute_ftilde`
- [GRHayL/Reconstruction/PPM/shock_detection_ftilde.c](../../../GRHayL/Reconstruction/PPM/shock_detection_ftilde.c) - `ghl_shock_detection_ftilde`
- [GRHayL/Reconstruction/PPM/slope_limit.c](../../../GRHayL/Reconstruction/PPM/slope_limit.c) - `ghl_slope_limit`
- [GRHayL/Reconstruction/PPM/steepen_var.c](../../../GRHayL/Reconstruction/PPM/steepen_var.c) - `ghl_steepen_var`
- [Unit_Tests/unit_test_ET_Legacy_reconstruction.c](../../../Unit_Tests/unit_test_ET_Legacy_reconstruction.c) - ET_Legacy PPM replay test

## See Also

- Depends on: [Face And Stencil Conventions](face-and-stencil-conventions.md)
- Depends on: [Struct Pack And Unpack](../../core/struct-pack-unpack.md)
- Validated by: [Reconstruction Tests And Reference Data](tests-and-reference-data.md)
- See also: [GRHayLib Runtime Parameters](../../implementations/grhaylib-runtime-parameters.md)
