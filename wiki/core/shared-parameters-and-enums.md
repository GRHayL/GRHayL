# Shared Parameters And Enums

Purpose: route GRHayL shared public parameters, enums, and cross-gem structs to
their owning source and KB pages. This page is a map, not a replacement for
public headers, Core packers, EOS setup, Con2Prim recovery, Reconstruction, or
Induction docs.

Ground truth starts in
[`GRHayL/include/ghl.h`](../../GRHayL/include/ghl.h),
[`GRHayL/include/ghl_io.h`](../../GRHayL/include/ghl_io.h),
[`GRHayL/include/ghl_debug.h`](../../GRHayL/include/ghl_debug.h), and
[`GRHayL/include/ghl_metric_helpers.h`](../../GRHayL/include/ghl_metric_helpers.h).
Use this page instead of creating separate public-header-install or
extrinsic-curvature-only pages.

## Public Header Surface

- [`GRHayL/include/ghl.h`](../../GRHayL/include/ghl.h) declares shared Core
  structs, public enums, EOS setup functions, Core packing helpers, metric
  helpers, stress-energy helpers, and small utility declarations.
- [`GRHayL/include/ghl_io.h`](../../GRHayL/include/ghl_io.h) declares
  `ghl_info`, `ghl_Warn_Error`, `ghl_error_keys`, and the terminating or
  warning IO macros.
- [`GRHayL/include/ghl_debug.h`](../../GRHayL/include/ghl_debug.h) provides
  inline primitive/conservative debug printers and wrapper macros.
- [`GRHayL/include/ghl_metric_helpers.h`](../../GRHayL/include/ghl_metric_helpers.h)
  provides inline raise/lower and vector-square helpers for 3D and 4D metric
  arrays.

## Enum Routes

| Public enum | Primary route | Use for |
| --- | --- | --- |
| `ghl_error_codes_t` | [errors, IO, debug, utilities](errors-io-debug-utilities.md) | Return-code vocabulary for EOS, Con2Prim, table, HDF5, atmosphere/floor, velocity, allocation, and utility errors. Message text is owned by [`abort_if_error.c`](../../GRHayL/GRHayL_Core/abort_if_error.c). |
| `ghl_con2prim_id_t` | [Con2Prim solver matrix](../gems/con2prim/solver-matrix.md) | Solver IDs such as `Noble2D`, `Font1D`, `Palenzuela1D`, `Newman1D`, entropy variants, and `None` backup sentinel. |
| `ghl_eos_t` | [EOS initialization and dispatch](../gems/eos/initialization-and-dispatch.md) | EOS family selection: simple, hybrid, or tabulated. |
| `ghl_eos_table_t` | [EOS initialization and dispatch](../gems/eos/initialization-and-dispatch.md) and tabulated EOS pages | Table-family selection, currently including `ghl_eos_table_stellarcollapse` and unknown/type-count sentinels. |
| `ghl_error_keys` | [errors, IO, debug, utilities](errors-io-debug-utilities.md) | IO macro control values in `ghl_io.h`, separate from recoverable `ghl_error_codes_t` returns. |

## `ghl_parameters`

`ghl_parameters` lives in
[`GRHayL/include/ghl.h`](../../GRHayL/include/ghl.h). Core initializes it in
[`GRHayL/GRHayL_Core/initialize_params.c`](../../GRHayL/GRHayL_Core/initialize_params.c).

Route fields by ownership:

- `main_routine` and `backup_routine[3]`: route `ghl_con2prim_id_t` values to
  [Con2Prim solver matrix](../gems/con2prim/solver-matrix.md) and recovery
  order to [Con2Prim recovery flow](../gems/con2prim/recovery-flow.md).
- `evolve_entropy`, `evolve_temp`, and `calc_prim_guess`: route recovery effects
  to [Con2Prim recovery flow](../gems/con2prim/recovery-flow.md) and EOS field
  coupling to [EOS initialization and dispatch](../gems/eos/initialization-and-dispatch.md).
- `max_Lorentz_factor` and `inv_sq_max_Lorentz_factor`: route velocity limiting
  to [Core velocity/u0 contract](velocity-u0-contract.md) and primitive-limit
  wrapper behavior to [Con2Prim limits and conversions](../gems/con2prim/limits-and-conversions.md).
- `psi6threshold`: route conservative/floor adjustment behavior to
  [Con2Prim limits and conversions](../gems/con2prim/limits-and-conversions.md).
- `con2prim_max_iterations` and `con2prim_solver_tolerance`: route solver
  behavior to [Con2Prim recovery flow](../gems/con2prim/recovery-flow.md) and
  method support to [Con2Prim solver matrix](../gems/con2prim/solver-matrix.md).
- `ppm_flattening_epsilon`, `ppm_flattening_omega1`,
  `ppm_flattening_omega2`, `ppm_shock_k0`, `ppm_shock_eta1`,
  `ppm_shock_eta2`, and `ppm_shock_epsilon`: route to
  [PPM flow](../gems/reconstruction/ppm-flow.md).
- `Lorenz_damping_factor`: route scalar-potential damping to
  [Induction gauge RHS contract](../gems/induction/gauge-rhs-contract.md).

`ghl_initialize_params` copies caller-selected control values, derives
`inv_sq_max_Lorentz_factor`, and sets default Con2Prim and `ppm_` controls.
Default values belong in
[`initialize_params.c`](../../GRHayL/GRHayL_Core/initialize_params.c) and the
existing [struct pack/unpack contract](struct-pack-unpack-contract.md), not in
new duplicated tables.

## `ghl_eos_parameters`

`ghl_eos_parameters` is declared in
[`GRHayL/include/ghl.h`](../../GRHayL/include/ghl.h). EOS initialization owns
its population; route setup details to
[EOS initialization and dispatch](../gems/eos/initialization-and-dispatch.md).

Use field groups as routes:

- EOS selection: `eos_type`, `table_type`, and `clean_sound_speed` route to EOS
  initialization and tabulated EOS pages.
- Atmosphere and floors: `rho_atm`, `rho_min`, `rho_max`, `tau_atm`,
  `press_atm`, `press_min`, `press_max`, `Y_e_atm`, `Y_e_min`, `Y_e_max`,
  `T_atm`, `T_min`, `T_max`, `eps_atm`, `eps_min`, `eps_max`,
  `entropy_atm`, `entropy_min`, and `entropy_max` route to
  [Atmosphere prescription contract](../gems/atmosphere/prescription-contract.md)
  and [Con2Prim limits and conversions](../gems/con2prim/limits-and-conversions.md).
- Hybrid/piecewise-polytrope fields, including `neos`, `rho_ppoly`,
  `Gamma_ppoly`, `K_ppoly`, `eps_integ_const`, `p_ppoly`, and `Gamma_th`,
  route to [hybrid piecewise-polytrope EOS](../gems/eos/hybrid-piecewise-polytrope.md).
- Tabulated table dimensions, arrays, bounds, interpolation auxiliaries,
  beta-equilibrium arrays, and `root_finding_precision` route to
  [tabulated interpolation and bounds](../gems/eos/tabulated-interpolation-and-bounds.md)
  and [tabulated EOS table contract](../gems/eos/tabulated-table-contract.md).

Core EOS setup entry points are declared in
[`GRHayL/include/ghl.h`](../../GRHayL/include/ghl.h) and implemented in
[`GRHayL/GRHayL_Core/initialize_eos.c`](../../GRHayL/GRHayL_Core/initialize_eos.c).

## Shared Struct Routes

| Struct | Route | Primary initializer |
| --- | --- | --- |
| `ghl_primitive_quantities` | [struct pack/unpack contract](struct-pack-unpack-contract.md), [Core velocity/u0 contract](velocity-u0-contract.md), [physics variables](../physics/variables-and-conventions.md) | [`initialize_primitives.c`](../../GRHayL/GRHayL_Core/initialize_primitives.c) |
| `ghl_conservative_quantities` | [struct pack/unpack contract](struct-pack-unpack-contract.md), [Con2Prim recovery flow](../gems/con2prim/recovery-flow.md), [physics variables](../physics/variables-and-conventions.md) | [`initialize_conservatives.c`](../../GRHayL/GRHayL_Core/initialize_conservatives.c) |
| `ghl_metric_quantities` | [metric/ADM contract](metric-adm-contract.md) | [`initialize_metric.c`](../../GRHayL/GRHayL_Core/initialize_metric.c) |
| `ghl_ADM_aux_quantities` | [metric/ADM contract](metric-adm-contract.md) | [`compute_ADM_auxiliaries.c`](../../GRHayL/GRHayL_Core/compute_ADM_auxiliaries.c) |
| `ghl_extrinsic_curvature` | [struct pack/unpack contract](struct-pack-unpack-contract.md) and source-term routes that need `K_{ij}` | [`initialize_extrinsic_curvature.c`](../../GRHayL/GRHayL_Core/initialize_extrinsic_curvature.c) |
| `ghl_stress_energy` | [stress-energy/smallb contract](stress-energy-smallb-contract.md) | [`initialize_stress_energy.c`](../../GRHayL/GRHayL_Core/initialize_stress_energy.c) |

These structs are public ABI surface. Do not duplicate full struct definitions
in KB pages; link
[`GRHayL/include/ghl.h`](../../GRHayL/include/ghl.h) for layout and use owner
pages above for behavioral contracts.

## Core Packers And Helpers

- [`GRHayL/GRHayL_Core/initialize_params.c`](../../GRHayL/GRHayL_Core/initialize_params.c)
  owns `ghl_parameters` initialization and default Con2Prim/PPM controls.
- [`GRHayL/GRHayL_Core/initialize_primitives.c`](../../GRHayL/GRHayL_Core/initialize_primitives.c)
  copies primitive fields; `u0` is not initialized there.
- [`GRHayL/GRHayL_Core/initialize_conservatives.c`](../../GRHayL/GRHayL_Core/initialize_conservatives.c)
  copies conservative fields using caller storage convention.
- [`GRHayL/GRHayL_Core/initialize_metric.c`](../../GRHayL/GRHayL_Core/initialize_metric.c)
  copies lapse, shift, and symmetric 3-metric inputs, then derives inverse and
  determinant fields.
- [`GRHayL/GRHayL_Core/initialize_extrinsic_curvature.c`](../../GRHayL/GRHayL_Core/initialize_extrinsic_curvature.c)
  copies six symmetric `K_{ij}` inputs into `ghl_extrinsic_curvature`.
- [`GRHayL/GRHayL_Core/initialize_stress_energy.c`](../../GRHayL/GRHayL_Core/initialize_stress_energy.c)
  copies ten symmetric tensor components into `ghl_stress_energy`.

## Review Hazards

- Con2Prim IDs are shared public API: adding or renaming an enum value requires
  matching solver support, name mapping, dispatch, build-list, and test review.
- EOS fields are shared across EOS, Atmosphere, Con2Prim, Flux_Source, and
  Neutrinos. Do not treat `ghl_eos_parameters` as Core-only.
- `ppm_` fields in `ghl_parameters` affect Reconstruction PPM behavior even
  though initialization lives in Core.
- `Lorenz_damping_factor` is stored in Core parameters but consumed by
  Induction gauge RHS callers.
- Error-code additions require both enum and `ghl_abort_if_error` message-map
  review.

## Repo-Local References

- [`GRHayL/include/ghl.h`](../../GRHayL/include/ghl.h)
- [`GRHayL/include/ghl_io.h`](../../GRHayL/include/ghl_io.h)
- [`GRHayL/include/ghl_debug.h`](../../GRHayL/include/ghl_debug.h)
- [`GRHayL/include/ghl_metric_helpers.h`](../../GRHayL/include/ghl_metric_helpers.h)
- [`GRHayL/GRHayL_Core/initialize_params.c`](../../GRHayL/GRHayL_Core/initialize_params.c)
- [`GRHayL/GRHayL_Core/initialize_primitives.c`](../../GRHayL/GRHayL_Core/initialize_primitives.c)
- [`GRHayL/GRHayL_Core/initialize_conservatives.c`](../../GRHayL/GRHayL_Core/initialize_conservatives.c)
- [`GRHayL/GRHayL_Core/initialize_metric.c`](../../GRHayL/GRHayL_Core/initialize_metric.c)
- [`GRHayL/GRHayL_Core/initialize_extrinsic_curvature.c`](../../GRHayL/GRHayL_Core/initialize_extrinsic_curvature.c)
- [`GRHayL/GRHayL_Core/initialize_stress_energy.c`](../../GRHayL/GRHayL_Core/initialize_stress_energy.c)
