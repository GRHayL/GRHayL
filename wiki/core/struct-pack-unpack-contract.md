# Struct Pack/Unpack Contract

This page routes Core struct and enum packing/unpacking questions. Ground truth
is `GRHayL/include/ghl.h`, the Core source files named below, and
`docs/raw/GRHayL_Core.dox`.

## Public Data Routes

- Core error, Con2Prim, EOS, and EOS-table enums live in
  `GRHayL/include/ghl.h`: `ghl_error_codes_t`, `ghl_con2prim_id_t`,
  `ghl_eos_t`, and `ghl_eos_table_t`. Con2Prim routine ownership is shared;
  use `wiki/public-api-map.md` before following non-Core implementations.
- `ghl_parameters` routes here for Core parameter packing, with Con2Prim and
  PPM defaults set by `GRHayL/GRHayL_Core/initialize_params.c`.
- `ghl_primitive_quantities` routes here for field packing/unpacking, except
  `u0`; `GRHayL/GRHayL_Core/initialize_primitives.c` leaves `u0` untouched and
  `GRHayL/GRHayL_Core/return_primitives.c` has no `u0` output.
- `ghl_conservative_quantities` routes here for field packing/unpacking in
  `GRHayL/GRHayL_Core/initialize_conservatives.c` and
  `GRHayL/GRHayL_Core/return_conservatives.c`.
- `ghl_extrinsic_curvature` routes here for symmetric field packing in
  `GRHayL/GRHayL_Core/initialize_extrinsic_curvature.c`.
- `ghl_metric_quantities` and `ghl_ADM_aux_quantities` route to
  `wiki/core/metric-adm-contract.md`; see `GRHayL/include/ghl.h`.
- `ghl_stress_energy` routes to
  `wiki/core/stress-energy-smallb-contract.md`; see `GRHayL/include/ghl.h`.
- `ghl_eos_parameters` routes to EOS initialization pages; its public layout is
  declared in `GRHayL/include/ghl.h`, but Core EOS dispatch is outside this
  page's ownership.

## Field-Copy Packers

Use these as copy contracts, not physics transforms:

- `ghl_initialize_primitives` copies density, pressure, epsilon, velocity,
  magnetic field, entropy, electron fraction, and temperature into
  `ghl_primitive_quantities`; it does not write `u0`
  (`GRHayL/GRHayL_Core/initialize_primitives.c`).
- `ghl_initialize_conservatives` copies density, energy, momentum, entropy, and
  electron fraction into `ghl_conservative_quantities`
  (`GRHayL/GRHayL_Core/initialize_conservatives.c`).
- `ghl_return_primitives` copies primitive fields back to caller pointers and
  does not return `u0` (`GRHayL/GRHayL_Core/return_primitives.c`).
- `ghl_return_conservatives` copies conservative fields back to caller pointers
  (`GRHayL/GRHayL_Core/return_conservatives.c`).
- `ghl_initialize_extrinsic_curvature` copies six inputs into symmetric
  `K[i][j]` slots (`GRHayL/GRHayL_Core/initialize_extrinsic_curvature.c`).

All functions in this section return `void`, dereference caller-provided
pointers, and perform no null, alias, range, or finiteness checks. They copy
values without changing densitization or physical representation. The
conservative struct is explicitly representation-neutral: its fields may hold
densitized or undensitized values, and the caller must track which convention
is present (`GRHayL/GRHayL_Core/initialize_conservatives.c` and
`GRHayL/GRHayL_Core/return_conservatives.c`).

Current `ghl_return_primitives` writes `*epsilon = prims->eps` twice. Both
writes use the same value, so this is not a second output or transform; retain
it as a source-review anomaly rather than a distinct API contract
(`GRHayL/GRHayL_Core/return_primitives.c`).

## Derived-Data Packers

- `ghl_initialize_params` copies caller-selected routines and booleans, computes
  `inv_sq_max_Lorentz_factor`, and sets default Con2Prim and PPM controls
  (`GRHayL/GRHayL_Core/initialize_params.c`).
- Its Con2Prim defaults are `con2prim_max_iterations = 30` and
  `con2prim_solver_tolerance = 1e-8`
  (`GRHayL/GRHayL_Core/initialize_params.c`).
- Its PPM defaults are `ppm_flattening_epsilon = 0.33`,
  `ppm_flattening_omega1 = 0.75`, `ppm_flattening_omega2 = 10.0`,
  `ppm_shock_k0 = 0.1`, `ppm_shock_eta1 = 20.0`,
  `ppm_shock_eta2 = 0.05`, and `ppm_shock_epsilon = 0.01`
  (`GRHayL/GRHayL_Core/initialize_params.c`).
- Metric initialization computes derived metric data; route
  `ghl_initialize_metric`, `ghl_compute_ADM_auxiliaries`, and
  `ghl_enforce_detgtij_and_initialize_ADM_metric` to
  `wiki/core/metric-adm-contract.md`.
- Stress-energy initialization and compute/return contracts route to
  `wiki/core/stress-energy-smallb-contract.md`.

`ghl_initialize_params` returns `void`, reads exactly three backup-routine
entries, computes `inv_sq_max_Lorentz_factor` by direct division, and performs
no validation. A null/short backup array, null output, zero or non-finite
Lorentz-factor input, or unsupported enum is outside its checked contract
(`GRHayL/include/ghl.h`; `GRHayL/GRHayL_Core/initialize_params.c`).

## Test Routes

- Core suite coverage around Core initialization and public helpers starts in
  `Unit_Tests/unit_test_grhayl_core_test_suite.c`.
- Reference-data generation for Core metric fixtures starts in
  `Unit_Tests/data_gen/unit_test_data_grhayl_core_test_suite.c`.
- The direct Core suite does not assert `ghl_initialize_params` defaults or a
  complete pack/return round trip. Legacy and generated-data paths exercise
  several helpers indirectly; use [Core tests and fixtures](tests-and-fixtures.md)
  for bounded coverage status.
