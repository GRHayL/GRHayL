# Velocity And `u0`

> Core `ghl_limit_v_and_compute_u0` velocity-cap rescaling, inward rounding retries, `u0` accuracy limits, and error returns. · Status: confirmed
> Up: [GRHayL Core](index.md)

## Summary

`ghl_limit_v_and_compute_u0` enforces the maximum Lorentz-factor speed cap from `params->inv_sq_max_Lorentz_factor` on transport velocities `prims->vU[0..2]` and computes `prims->u0`, returning `ghl_error_u0_singular` when a required cap is non-subluminal, when the final Eulerian speed squared is nonfinite, luminal/superluminal, or above the cap, or when `u0` is nonfinite. EOS-dependent primitive floors belong separately to Con2Prim's `ghl_enforce_primitive_limits_and_compute_u0`.

## Detail

### Inputs, Rescaling Retries, And Floating-Point Accuracy Limits

`ghl_limit_v_and_compute_u0` (`limit_v_and_compute_u0.c`) expects initialized `ghl_parameters` (finite `max_Lorentz_factor >= 1`), initialized ADM `ghl_metric_quantities` (finite positive lapse, shift `betaU`, positive-definite `gammaDD`), transport velocities `prims->vU[0..2]`, and caller-initialized `*speed_limited` (which it OR-updates on limiting and never resets to `false`).

- It computes the Eulerian squared speed via `ghl_compute_vec2_from_vec3D`. An already subluminal state below the cap keeps its velocity unchanged regardless of extreme cap settings.
- When limiting is required, it returns `ghl_error_u0_singular` **before** mutating `vU` if the squared-speed bound is nonfinite or not strictly subluminal (`>= 1`).
- After rescaling `vU[0..2]`, it recomputes the achieved metric speed. If double-precision rounding leaves the achieved speed above the cap, it retries inward rescaling up to 64 times with a doubling inward margin and strictly inward-rounded scale (spanning the bound within 54 doublings; at very small lapse near `v^i = -beta^i`, quantization can place the accepted speed well inside the cap, down to `u0 = 1/alpha`).
- It returns `ghl_error_u0_singular` if the final achieved speed is nonfinite, luminal/superluminal, or above the cap, or if computed `u0` is nonfinite; otherwise it writes `prims->u0` and returns `ghl_success`.
- **Double-precision accuracy:** for limited states, relative `u0` accuracy is roughly `DBL_EPSILON * max_Lorentz_factor^2` (~`2e-6` at `1e5`, ~`2e-2` at `1e7`). Near `5e7`, a double-precision accepted state can be superluminal in exact arithmetic; above ~`1.34e8`, `1 - inv_sq_max_Lorentz_factor` rounds to `1.0`, so any state requiring limiting returns `ghl_error_u0_singular`. Pointers, metric signature, and lapse positivity are not checked.

Evidence (descriptive behavior): `limit_v_and_compute_u0.c::ghl_limit_v_and_compute_u0`; corroboration: `ghl.h::ghl_error_u0_singular`.

### Separation From Con2Prim Limiting And Test Coverage

`ghl_enforce_primitive_limits_and_compute_u0` (`GRHayL/Con2Prim/enforce_primitive_limits_and_compute_u0.c`) applies EOS-specific primitive floors/ceilings first and then calls Core `ghl_limit_v_and_compute_u0` as its final step; see [Limits And Conversions](../gems/con2prim/limits-and-conversions.md).

- `Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c` tests the Con2Prim wrapper against reference data and calls `ghl_limit_v_and_compute_u0` directly for ordinary-cap success (including shifted `W_max = 10` rounding), extreme rounded-cap rejection, achieved-speed enforcement, and subluminal no-limit passthrough.
- `Unit_Tests/unit_test_code_error.c` key `4` calls `ghl_limit_v_and_compute_u0` with `NaN` velocities to verify `ghl_error_u0_singular`.
- `unit_test_hybrid_flux.c` and `unit_test_tabulated_flux.c` call `ghl_limit_v_and_compute_u0` on sub-cap face states only.

Evidence (descriptive behavior): `unit_test_enforce_primitive_limits_and_compute_u0.c::ghl_limit_v_and_compute_u0`; corroboration: `unit_test_code_error.c::ghl_error_u0_singular`.

## Sources

- [GRHayL/GRHayL_Core/limit_v_and_compute_u0.c](../../GRHayL/GRHayL_Core/limit_v_and_compute_u0.c) - velocity cap, inward retry loop, and `u0` calculation
- [GRHayL/include/ghl.h](../../GRHayL/include/ghl.h) - `ghl_limit_v_and_compute_u0`, `ghl_parameters`, `ghl_error_u0_singular`
- [GRHayL/include/ghl_metric_helpers.h](../../GRHayL/include/ghl_metric_helpers.h) - `ghl_compute_vec2_from_vec3D`
- [Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c](../../Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c) - direct Core limiter and wrapper checks
- [Unit_Tests/unit_test_code_error.c](../../Unit_Tests/unit_test_code_error.c) - key `4` singular-`u0` test

## See Also

- Depends on: [Metric And ADM Helpers](metric-adm.md)
- Depends on: [Shared Parameters And Enums](shared-parameters-and-enums.md)
- See also: [Limits And Conversions](../gems/con2prim/limits-and-conversions.md)
- Validated by: [Core Tests And Reference Data](tests-and-reference-data.md)
