# Face And Stencil Conventions

> Left-face `Ur`/`Ul` orientation, 4-/5-/6-point stencil layouts, centered-helper mapping, and caller assumptions across `GRHayL/Reconstruction/`. · Status: confirmed
> Up: [Reconstruction](index.md)

## Summary

All standard public PLM, PPM, and WENO-z reconstruction wrappers return the right (`Ur`) and left (`Ul`) reconstructed states on the **left face `i-1/2`** of cell `i`, using uniform-grid stencils (`U[4]` for PLM; `U[6]` for PPM and WENO-z). By contrast, the 5-point PPM and WENO-z helper routines compute both faces around a single stencil center `c` (`Ur` at `c+1/2` left side, `Ul` at `c-1/2` right side). `docs/raw/Reconstruction.dox` stays read-only during KB updates.

## Detail

### Left-Face Wrappers Versus Centered 5-Point Helpers

For standard public wrappers (`ghl_minmod_reconstruction`, `ghl_mc_reconstruction`, `ghl_superbee_reconstruction`, `ghl_ppm_reconstruction`, `ghl_ppm_reconstruction_with_steepening`, `ghl_wenoz_reconstruction`), face `i-1/2` uses:
- `Ur = U(i - 1/2 + epsilon)` (right side of face `i-1/2`, current-cell side);
- `Ul = U(i - 1/2 - epsilon)` (left side of face `i-1/2`, neighbor-cell side).

| Method / layer | Input width | Cell span | Notes |
| --- | --- | --- | --- |
| PLM public wrappers | `U[4]` | `i-2 .. i+1` | Tests pass `&var[index - 2]`. (Ev: `GRHayL/include/ghl_reconstruction.h::ghl_minmod_reconstruction`) |
| PPM & WENO-z public wrappers | `U[6]` | `i-3 .. i+2` | Tests pass `&var[index - 3]`; `ghl_compute_ftilde` uses aligned `pressure[6]` and `v_flux_dirn[6]`. (Ev: `GRHayL/include/ghl_reconstruction.h::ghl_wenoz_reconstruction`) |
| PPM & WENO-z centered helpers | `U[5]` | `c-2 .. c+2` (`MINUS2..PLUS_2`) | Writes left state at `c+1/2` to `Ur` and right state at `c-1/2` to `Ul`; wrapper calls on `c = i-1` (taking `Ur -> Ul`) and `c = i` (taking `Ul -> Ur`). (Ev: `GRHayL/include/ghl_reconstruction.h::reconstruction_stencil`) |

Evidence (scientific/interface requirement): `GRHayL/include/ghl_reconstruction.h::reconstruction_stencil`; corroboration: `GRHayL/Reconstruction/WENOZ/wenoz_reconstruction.c` and `GRHayL/Reconstruction/PPM/ppm_reconstruction.c`.

### Exposure, Uniform-Grid Assumptions, And Caller Impact

- **Header exposure:** `GRHayL/include/ghl_reconstruction.h` (installed in default and `--disable-hdf5` builds) declares both standard wrappers and internal-support helpers (`ghl_minmod`, `ghl_maxmod`, `ghl_shock_detection_ftilde`, `ghl_slope_limit`, `ghl_steepen_var`, plus 5-point PPM/WENO-z helpers) with global linkage; only WENO-z's static `mc` helper is translation-unit private.
- **No production in-repo callers:** Only `Unit_Tests/` and `Unit_Tests/data_gen/` call Reconstruction directly; `Flux_Source`, `Induction`, and `GRHayLib` contain no in-tree call edges (downstream host codes connect reconstructed states to fluxes).
- **Numerical assumptions:** Fixed coefficients assume uniform coordinate spacing. Routines take no grid widths, metric, or EOS struct, return no error status, and do not enforce positivity, EOS consistency, or relativistic velocity bounds. Changing PPM parameters in `ghl_parameters` (`ppm_flattening_*`, `ppm_shock_*`) affects `ghl_compute_ftilde`, `ghl_shock_detection_ftilde`, `ghl_steepen_var`, ET_Legacy reconstruction, and `GRHayLib/param.ccl`.

Evidence (descriptive behavior): `GRHayL/Reconstruction/make.code.defn`; corroboration: `docs/raw/Reconstruction.dox`.

## Sources

- [GRHayL/include/ghl_reconstruction.h](../../../GRHayL/include/ghl_reconstruction.h) - wrapper/helper declarations and `enum reconstruction_stencil`
- [GRHayL/Reconstruction/make.code.defn](../../../GRHayL/Reconstruction/make.code.defn) - `PLM`, `PPM`, `WENOZ` subdirectories
- [GRHayL/Reconstruction/PPM/ppm_reconstruction.c](../../../GRHayL/Reconstruction/PPM/ppm_reconstruction.c) - 6-point to 5-point wrapper mapping
- [GRHayL/Reconstruction/WENOZ/wenoz_reconstruction.c](../../../GRHayL/Reconstruction/WENOZ/wenoz_reconstruction.c) - 6-point to 5-point WENO-z mapping
- [docs/raw/Reconstruction.dox](../../../docs/raw/Reconstruction.dox) - read-only Reconstruction documentation

## See Also

- Implements: [PLM Limiters](plm-limiters.md)
- Implements: [PPM Flow](ppm-flow.md)
- Implements: [WENO-z Method](wenoz-method.md)
- Validated by: [Reconstruction Tests And Reference Data](tests-and-reference-data.md)
