# WENO-z Method

> WENO-z 6-point left-face wrapper, 5-point centered helper, fixed weights and MC blend, numerical overflow/positivity caveats, and CI coverage. · Status: confirmed
> Up: [Reconstruction](index.md)

## Summary

`GRHayL/Reconstruction/WENOZ/` provides the 6-point left-face wrapper `ghl_wenoz_reconstruction` backed by two shifted calls to the 5-point centered helper `ghl_wenoz_reconstruction_right_left_faces`. Weights and regularization are fixed in source, extreme smoothness contrasts can overflow weights, positivity is not guaranteed even for strictly positive stencils, and WENO-z replay runs in CI workflows rather than `.github/run_tests.sh`.

## Detail

### Stencil Adaptation, Fixed Constants, And Numerical Limits

- **Wrapper-to-helper mapping:** `ghl_wenoz_reconstruction` (`wenoz_reconstruction.c`) takes `U[6]` (`i-3 .. i+2`), calls `ghl_wenoz_reconstruction_right_left_faces` on `&U[0]` (assigning helper right-face output to wrapper `Ul`), and calls it on `&U[1]` (assigning helper left-face output to wrapper `Ur`).
- **Fixed parameters:** No runtime WENO parameter struct exists. `wenoz_reconstruction_right_left_faces.c` (which also holds upstream provenance/license comments) hard-codes linear weights `{0.1, 0.6, 0.3}`, regularization epsilon `1e-100`, and MC fallback coefficient `2.0`.
- **Overflow and non-positivity caveats:** Fixed regularization with unscaled double-precision weight ratios can overflow on extreme smoothness contrasts (for example, `{0, 0, 0, 1e105, 1e105}` yields nonfinite outputs). The MC blend does not guarantee positivity or local-bound preservation: on the strictly positive centered stencil `{9, 1, 1, 9, 25}`, equal smoothness indicators and zero MC slope retain a WENO left-face value near `-1/3`. Callers enforce physical floors separately.

Evidence (descriptive behavior): `GRHayL/Reconstruction/WENOZ/wenoz_reconstruction_right_left_faces.c::ghl_wenoz_reconstruction_right_left_faces`; corroboration: `GRHayL/Reconstruction/WENOZ/wenoz_reconstruction.c::ghl_wenoz_reconstruction`.

### Build Membership And Test Routing

`GRHayL/Reconstruction/WENOZ/make.code.defn` builds both WENO-z C files. `Unit_Tests/unit_test_WENOZ_reconstruction.c` compares wrapper `Ur`/`Ul` against `WENOZ_reconstruction_{input,output,output_pert}.bin` (produced by `unit_test_data_WENOZ_reconstruction.c` over `3 <= index < arraylength - 2` with zero boundary placeholders) and does not call the 5-point helper directly. Every compiler workflow runs `WENOZ_reconstruction` in its `reconstruction` job matrix, whereas `.github/run_tests.sh` omits it.

Evidence (descriptive behavior): `GRHayL/Reconstruction/WENOZ/make.code.defn`; corroboration: `Unit_Tests/unit_test_WENOZ_reconstruction.c`.

## Sources

- [GRHayL/include/ghl_reconstruction.h](../../../GRHayL/include/ghl_reconstruction.h) - WENO-z declarations
- [GRHayL/Reconstruction/WENOZ/make.code.defn](../../../GRHayL/Reconstruction/WENOZ/make.code.defn) - WENO-z build list
- [GRHayL/Reconstruction/WENOZ/wenoz_reconstruction.c](../../../GRHayL/Reconstruction/WENOZ/wenoz_reconstruction.c) - 6-point wrapper
- [GRHayL/Reconstruction/WENOZ/wenoz_reconstruction_right_left_faces.c](../../../GRHayL/Reconstruction/WENOZ/wenoz_reconstruction_right_left_faces.c) - 5-point helper, weights, and MC blend
- [Unit_Tests/unit_test_WENOZ_reconstruction.c](../../../Unit_Tests/unit_test_WENOZ_reconstruction.c) - direct WENO-z replay test
- [Unit_Tests/data_gen/unit_test_data_WENOZ_reconstruction.c](../../../Unit_Tests/data_gen/unit_test_data_WENOZ_reconstruction.c) - WENO-z reference generator

## See Also

- Depends on: [Face And Stencil Conventions](face-and-stencil-conventions.md)
- Validated by: [Reconstruction Tests And Reference Data](tests-and-reference-data.md)
- Contrasts with: [PPM Flow](ppm-flow.md)
