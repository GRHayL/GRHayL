# Struct Pack And Unpack

> Core field-copy packers (`ghl_initialize_*`, `ghl_return_*`), `ghl_initialize_params` defaults, and unchecked preconditions. · Status: confirmed
> Up: [GRHayL Core](index.md)

## Summary

Core packing and unpacking routines copy primitive, conservative, and extrinsic-curvature scalars into and out of `ghl.h` structs without changing physical normalization or densitization, while `ghl_initialize_params` computes `inv_sq_max_Lorentz_factor` and sets default Con2Prim and PPM controls. Every routine here returns `void` and performs no pointer, range, or finiteness checks.

## Detail

### Field-Copy Packers And Unpackers

| Routine | File | Behavior & notable quirks |
| --- | --- | --- |
| `ghl_initialize_primitives` | `initialize_primitives.c` | Copies `rho`, `press`, `eps`, `vU[0..2]`, `BU[0..2]`, `entropy`, `Y_e`, `temperature`; **leaves `u0` untouched**. (Ev: `initialize_primitives.c::ghl_initialize_primitives`) |
| `ghl_return_primitives` | `return_primitives.c` | Copies primitive fields to caller pointers; omits `u0`; currently assigns `*epsilon = prims->eps` twice identically (harmless duplicate assignment). (Ev: `return_primitives.c::ghl_return_primitives`) |
| `ghl_initialize_conservatives` / `ghl_return_conservatives` | `initialize_conservatives.c`, `return_conservatives.c` | Copies `rho`, `tau`, `SD[0..2]`, `entropy`, `Y_e` without converting between densitized and undensitized representations. (Ev: `initialize_conservatives.c::ghl_initialize_conservatives`) |
| `ghl_initialize_extrinsic_curvature` | `initialize_extrinsic_curvature.c` | Copies six symmetric inputs into `K[3][3]`. (Ev: `initialize_extrinsic_curvature.c::ghl_initialize_extrinsic_curvature`) |

None of these routines validates pointers, aliasing, or finiteness; callers track whether `ghl_conservative_quantities` holds densitized or undensitized values.

### `ghl_initialize_params` Defaults And Preconditions

`ghl_initialize_params` (`initialize_params.c`) copies caller-supplied routine IDs and flags, reads all three entries of `backup_routine`, computes `inv_sq_max_Lorentz_factor = 1.0 / (max_Lorentz_factor * max_Lorentz_factor)` without validation (caller must supply finite `max_Lorentz_factor >= 1`), and sets these defaults:

- **Con2Prim:** `con2prim_max_iterations = 30`, `con2prim_solver_tolerance = 1e-8`.
- **PPM:** `ppm_flattening_epsilon = 0.33`, `ppm_flattening_omega1 = 0.75`, `ppm_flattening_omega2 = 10.0`, `ppm_shock_k0 = 0.1`, `ppm_shock_eta1 = 20.0`, `ppm_shock_eta2 = 0.05`, `ppm_shock_epsilon = 0.01`.

The direct Core test suite does not assert every `ghl_initialize_params` default or full pack/unpack round trips; see [Core Tests And Reference Data](tests-and-reference-data.md).

Evidence (descriptive behavior): `GRHayL/GRHayL_Core/initialize_params.c::ghl_initialize_params`; corroboration: `GRHayL/include/ghl.h::ghl_parameters`.

## Sources

- [GRHayL/include/ghl.h](../../GRHayL/include/ghl.h) - struct and packer declarations
- [GRHayL/GRHayL_Core/initialize_params.c](../../GRHayL/GRHayL_Core/initialize_params.c) - parameter defaults and `inv_sq_max_Lorentz_factor`
- [GRHayL/GRHayL_Core/initialize_primitives.c](../../GRHayL/GRHayL_Core/initialize_primitives.c) - primitive packing (leaves `u0` untouched)
- [GRHayL/GRHayL_Core/return_primitives.c](../../GRHayL/GRHayL_Core/return_primitives.c) - primitive unpacking and duplicate `*epsilon` write
- [GRHayL/GRHayL_Core/initialize_conservatives.c](../../GRHayL/GRHayL_Core/initialize_conservatives.c) - conservative packing
- [GRHayL/GRHayL_Core/return_conservatives.c](../../GRHayL/GRHayL_Core/return_conservatives.c) - conservative unpacking
- [GRHayL/GRHayL_Core/initialize_extrinsic_curvature.c](../../GRHayL/GRHayL_Core/initialize_extrinsic_curvature.c) - symmetric `K_{ij}` packing
- [docs/raw/GRHayL_Core.dox](../../docs/raw/GRHayL_Core.dox) - Core API documentation

## See Also

- Depends on: [Shared Parameters And Enums](shared-parameters-and-enums.md)
- Validated by: [Core Tests And Reference Data](tests-and-reference-data.md)
- See also: [Velocity And `u0`](velocity-u0.md)
- See also: [PPM Flow](../gems/reconstruction/ppm-flow.md)
