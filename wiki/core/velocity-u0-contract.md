# Velocity/u0 Contract

## Purpose

Route Core velocity limiting and `u0` computation. This page separates
`ghl_limit_v_and_compute_u0` in Core from the Con2Prim primitive-limit wrapper;
source, headers, Doxygen, and tests remain authority
(`GRHayL/GRHayL_Core/limit_v_and_compute_u0.c`,
`GRHayL/include/ghl.h`, `docs/raw/GRHayL_Core.dox`).

## Core Routine

`ghl_limit_v_and_compute_u0` is declared as a Core public routine taking
`ghl_parameters`, initialized ADM metric quantities, mutable primitives, and a
caller-provided `speed_limited` pointer (`GRHayL/include/ghl.h:496`).
`ghl_parameters` carries `inv_sq_max_Lorentz_factor`, and
`ghl_metric_quantities` carries `lapseinv`, `lapseinv2`, `betaU`, and
`gammaDD`, all used by the implementation (`GRHayL/include/ghl.h:102`,
`GRHayL/include/ghl.h:154`,
`GRHayL/GRHayL_Core/limit_v_and_compute_u0.c:21`).

Inputs expected before call:
- Initialized parameters from the normal Core parameter path, because the speed
  cap is read from `params->inv_sq_max_Lorentz_factor`
  (`GRHayL/include/ghl.h:102`,
  `GRHayL/GRHayL_Core/limit_v_and_compute_u0.c:25`).
- Initialized ADM metric from the Core metric path, because the routine reads
  lapse, shift, and spatial metric fields (`GRHayL/include/ghl.h:154`,
  `GRHayL/GRHayL_Core/limit_v_and_compute_u0.c:21`).
- Primitive velocity fields `prims->vU[0..2]`; these are the mutable transport
  velocities stored in `ghl_primitive_quantities`
  (`GRHayL/include/ghl.h:131`,
  `GRHayL/GRHayL_Core/limit_v_and_compute_u0.c:21`).
- Caller-initialized `speed_limited`. The routine OR-updates the bool when it
  applies the speed cap; it does not reset the flag to `false`
  (`GRHayL/GRHayL_Core/limit_v_and_compute_u0.c:32`).

Outputs and mutation:
- If the velocity norm exceeds the Lorentz-factor cap, the routine rescales and
  writes all three `prims->vU` components
  (`GRHayL/GRHayL_Core/limit_v_and_compute_u0.c:26`).
- It writes `prims->u0` after the optional velocity cap
  (`GRHayL/GRHayL_Core/limit_v_and_compute_u0.c:39`).
- It returns `ghl_error_u0_singular` if the computed `u0` is `NaN`; otherwise
  it returns `ghl_success` (`GRHayL/include/ghl.h:47`,
  `GRHayL/GRHayL_Core/limit_v_and_compute_u0.c:41`).
- The vector norm uses `ghl_compute_vec2_from_vec3D`, whose helper contract is
  the metric-indexed 3-vector square in `GRHayL/include/ghl_metric_helpers.h:158`
  (`GRHayL/GRHayL_Core/limit_v_and_compute_u0.c:22`).

## Not Con2Prim Limits

`ghl_enforce_primitive_limits_and_compute_u0` is a Con2Prim routine, not this
Core routine. The Con2Prim wrapper applies EOS-specific primitive floors,
ceilings, pressure/internal-energy updates, and then calls the Core
`ghl_limit_v_and_compute_u0` as its final velocity/u0 step; keep behavioral
claims for EOS limiting routed to
`wiki/gems/con2prim/limits-and-conversions.md` and
`GRHayL/Con2Prim/enforce_primitive_limits_and_compute_u0.c:10`,
`GRHayL/Con2Prim/enforce_primitive_limits_and_compute_u0.c:17`,
`GRHayL/Con2Prim/enforce_primitive_limits_and_compute_u0.c:53`,
`GRHayL/Con2Prim/enforce_primitive_limits_and_compute_u0.c:76`, and
`GRHayL/Con2Prim/enforce_primitive_limits_and_compute_u0.c:93`.
The public Core declaration for `ghl_limit_v_and_compute_u0` is in
`GRHayL/include/ghl.h:496`; Con2Prim test coverage calls the wrapper by name in
`Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c:182`.

## Tests

Normal-path routing: `Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c`
initializes parameters, hybrid EOS, metric, ADM auxiliaries, and primitives,
then calls `ghl_enforce_primitive_limits_and_compute_u0` and checks `prims.u0`
against fixture data (`Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c:35`,
`Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c:164`,
`Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c:182`,
`Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c:201`). Treat this
as wrapper-path coverage for the Core limiter, not standalone full coverage of
all Core inputs.

Error-path routing: `Unit_Tests/unit_test_code_error.c` initializes metric data,
sets primitive velocities to `NaN`, calls `ghl_limit_v_and_compute_u0`
directly, and maps test key 4 to `ghl_error_u0_singular`
(`Unit_Tests/unit_test_code_error.c:122`,
`Unit_Tests/unit_test_code_error.c:141`,
`Unit_Tests/unit_test_code_error.c:146`,
`Unit_Tests/unit_test_code_error.c:649`). This covers the singular-`u0` error
route without implying broader standalone test coverage.
