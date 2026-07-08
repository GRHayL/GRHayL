# GRHayL Core Hub

This page is a router for GRHayL_Core work, not replacement documentation.
Use it to pick focused Core KB pages and repo-local ground truth before
reviewing source, Doxygen, tests, or downstream integrations.

If this page conflicts with source, headers, tests, CI, or Doxygen source,
trust the underlying repo files and update the KB route.

## Focused Core Pages

| Page | Route there for |
| --- | --- |
| [Struct pack/unpack contract](struct-pack-unpack-contract.md) | `ghl_parameters`, primitive/conservative structs, extrinsic curvature packing, and field-copy return helpers from `ghl.h`. |
| [Metric/ADM contract](metric-adm-contract.md) | ADM metric initialization, determinant-enforced setup, ADM auxiliaries, and inline vector metric helpers. |
| [Velocity/u0 contract](velocity-u0-contract.md) | Core `ghl_limit_v_and_compute_u0`, velocity limiting, `u0`, and singular-`u0` error routing. |
| [Stress-energy/smallb contract](stress-energy-smallb-contract.md) | `smallb`, `b2`, stress-energy packing, and `Tmunu` compute/return helpers. |
| [EOS dispatch contract](eos-dispatch-contract.md) | Core EOS initialization wrappers, EOS/flux function-pointer dispatch, and HDF5-disabled tabulated routes. |
| [Errors/IO/debug/utilities](errors-io-debug-utilities.md) | Core logging, abort-on-error, min/max/clamp helpers, and debug-print headers. |
| [Tests and fixtures](tests-and-fixtures.md) | Core unit suite, fixture generation, perturbation checks, and adjacent Con2Prim/stress-energy tests. |

## Source Inventory

| Ground truth | Use it for |
| --- | --- |
| [`GRHayL/GRHayL_Core/make.code.defn`](../../GRHayL/GRHayL_Core/make.code.defn) | Compiled Core source list and source-file ownership audit. |
| [`GRHayL/include/ghl.h`](../../GRHayL/include/ghl.h) | Shared Core public structs, enums, initializer declarations, metric/stress-energy declarations, EOS setup declarations, and small utility declarations. |
| [`GRHayL/include/ghl_metric_helpers.h`](../../GRHayL/include/ghl_metric_helpers.h) | Header-only 3D/4D raise/lower and vector-square helpers. |
| [`GRHayL/include/ghl_io.h`](../../GRHayL/include/ghl_io.h) | Logging/error macros and `ghl_Warn_Error` API. |
| [`GRHayL/include/ghl_debug.h`](../../GRHayL/include/ghl_debug.h) | Header-only primitive/conservative debug printing helpers. |
| [`GRHayL/include/ghl_eos_functions.h`](../../GRHayL/include/ghl_eos_functions.h) | External declarations for EOS and flux function pointers. |
| [`GRHayL/include/ghl_eos_functions_declaration.h`](../../GRHayL/include/ghl_eos_functions_declaration.h) | Storage definitions for EOS and flux function pointers. |
| [`GRHayL/include/make.code.defn`](../../GRHayL/include/make.code.defn) | Public header install/build availability list. |

Related routers: [KB index](../index.md), [catalog](../catalog.md),
[source map](../source-map.md), and [public API map](../public-api-map.md).

## Boundary Audit

- Core-implemented declarations in [`GRHayL/include/ghl.h`](../../GRHayL/include/ghl.h)
  route to the focused Core pages above. Use the source-file ownership table
  below before treating a declaration as a Core implementation.
- `ghl_get_con2prim_routine_name` is declared in
  [`GRHayL/include/ghl.h`](../../GRHayL/include/ghl.h) but implemented in
  [`GRHayL/Con2Prim/get_con2prim_routine_name.c`](../../GRHayL/Con2Prim/get_con2prim_routine_name.c).
  Route it as a shared public declaration exposed by `ghl.h`, not as a Core
  implementation.
- `ghl_compute_SU_Bsq_Ssq_BdotS` is declared in
  [`GRHayL/include/ghl.h`](../../GRHayL/include/ghl.h) but implemented in
  [`GRHayL/Con2Prim/compute_SU_Bsq_Ssq_BdotS.c`](../../GRHayL/Con2Prim/compute_SU_Bsq_Ssq_BdotS.c).
  Route implementation details to Con2Prim, especially
  [Con2Prim limits and conversions](../gems/con2prim/limits-and-conversions.md).
- [`GRHayL/include/make.code.defn`](../../GRHayL/include/make.code.defn)
  controls which public headers are installed or otherwise made available by
  the build. Header availability claims should link that file.

## Source-File Ownership

Every source file listed by
[`GRHayL/GRHayL_Core/make.code.defn`](../../GRHayL/GRHayL_Core/make.code.defn)
routes to exactly one focused Core page:

| Core file | Focused node |
| --- | --- |
| `abort_if_error.c` | [Errors/IO/debug/utilities](errors-io-debug-utilities.md) |
| `compute_ADM_auxiliaries.c` | [Metric/ADM contract](metric-adm-contract.md) |
| `compute_smallb_and_b2.c` | [Stress-energy/smallb contract](stress-energy-smallb-contract.md) |
| `compute_TDNmunu.c` | [Stress-energy/smallb contract](stress-energy-smallb-contract.md) |
| `compute_TUPmunu.c` | [Stress-energy/smallb contract](stress-energy-smallb-contract.md) |
| `enforce_detgij_and_initialize_ADM_metric.c` | [Metric/ADM contract](metric-adm-contract.md) |
| `info_warn_error.c` | [Errors/IO/debug/utilities](errors-io-debug-utilities.md) |
| `initialize_conservatives.c` | [Struct pack/unpack contract](struct-pack-unpack-contract.md) |
| `initialize_eos.c` | [EOS dispatch contract](eos-dispatch-contract.md) |
| `initialize_extrinsic_curvature.c` | [Struct pack/unpack contract](struct-pack-unpack-contract.md) |
| `initialize_metric.c` | [Metric/ADM contract](metric-adm-contract.md) |
| `initialize_params.c` | [Struct pack/unpack contract](struct-pack-unpack-contract.md) |
| `initialize_primitives.c` | [Struct pack/unpack contract](struct-pack-unpack-contract.md) |
| `initialize_stress_energy.c` | [Stress-energy/smallb contract](stress-energy-smallb-contract.md) |
| `limit_v_and_compute_u0.c` | [Velocity/u0 contract](velocity-u0-contract.md) |
| `min_max_clamp.c` | [Errors/IO/debug/utilities](errors-io-debug-utilities.md) |
| `return_conservatives.c` | [Struct pack/unpack contract](struct-pack-unpack-contract.md) |
| `return_primitives.c` | [Struct pack/unpack contract](struct-pack-unpack-contract.md) |
| `return_stress_energy.c` | [Stress-energy/smallb contract](stress-energy-smallb-contract.md) |
