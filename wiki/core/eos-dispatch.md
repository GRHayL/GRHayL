# EOS Dispatch

> Core EOS setup wrappers, process-wide function-pointer globals, initialization ordering, and `GHL_DISABLE_HDF5` dispatch behavior. · Status: confirmed
> Up: [GRHayL Core](index.md)

## Summary

`GRHayL/GRHayL_Core/initialize_eos.c` holds the storage for GRHayL's process-wide EOS and Con2Prim function-pointer globals and implements the public simple, hybrid, and tabulated setup wrappers. Direct parameter initializers depend on prior function-pointer dispatch setup; wrappers assign global pointers before validating parameters without rollback, while tabulated initialization constructs its struct transactionally and returns `ghl_error_used_disabled_hdf5` under `GHL_DISABLE_HDF5`.

## Detail

### Process-Wide Function-Pointer Dispatch And Ordering

EOS dispatch uses process-wide mutable function-pointer globals rather than members of `ghl_eos_parameters`:
- `ghl_con2prim_multi_method` is declared in `ghl_con2prim.h` and defined in `initialize_eos.c`.
- EOS pointers are declared `extern` in `ghl_eos_functions.h` and defined in `initialize_eos.c` by including the non-installed internal file `ghl_eos_functions_declaration.h`. Static storage starts zero-initialized (not callable).

`ghl_initialize_eos_functions(eos_type)` executes in three stages:
1. Calls `NRPyEOS_initialize_hybrid_functions` (assigning hybrid pointers and `ghl_compute_h_and_cs2`).
2. In HDF5 builds, calls `NRPyEOS_initialize_tabulated_functions` (assigning tabulated pointers including `ghl_tabulated_free_beq_quantities`, temporarily overwriting `ghl_compute_h_and_cs2`).
3. For a recognized `ghl_eos_t`, sets `ghl_con2prim_multi_method` and restores `ghl_compute_h_and_cs2` to the requested family (`NRPyEOS_hybrid_compute_enthalpy_and_cs2` for simple/hybrid; tabulated for `ghl_eos_tabulated`).

An unrecognized `ghl_eos_t` has no `else` error path: `ghl_initialize_eos_functions` returns `void`, leaves `ghl_con2prim_multi_method` unchanged, and leaves `ghl_compute_h_and_cs2` pointing to hybrid (no-HDF5) or tabulated (HDF5). The deprecated generic `ghl_calculate_HLLE_fluxes_dirn*` pointer globals remain zero-initialized and unassigned; call direct named Flux_Source variants instead.

Direct parameter initializers (`ghl_initialize_simple_eos`, `ghl_initialize_hybrid_eos`, `ghl_initialize_tabulated_eos`) call those global pointers internally, so callers must call `ghl_initialize_eos_functions` first or use a `*_functions_and_params` wrapper. Wrappers mutate global dispatch before parameter validation (no rollback on error): simple/hybrid leave the destination struct untouched on error; tabulated frees any candidate and publishes an empty tagged aggregate on failure (callers must free live tables before reinitialization). For family-specific parameter validation, rectangular extrema, `tau_atm`, and tabulated wrapper defaults, see [EOS Initialization And Dispatch](../gems/eos/initialization-and-dispatch.md).

Evidence (descriptive behavior): `GRHayL/GRHayL_Core/initialize_eos.c::ghl_initialize_eos_functions`; corroboration: `GRHayL/include/ghl_eos_functions_declaration.h`.

### `GHL_DISABLE_HDF5` Dispatch Divergence And Core Test Coverage

- **No-HDF5 divergence:** Under `GHL_DISABLE_HDF5`, `ghl_initialize_tabulated_eos` and `ghl_initialize_tabulated_eos_functions_and_params` return `ghl_error_used_disabled_hdf5` (publishing an empty `ghl_eos_tabulated` aggregate on valid non-live storage) without modifying dispatch. In contrast, calling `ghl_initialize_eos_functions(ghl_eos_tabulated)` directly hits `GHL_HDF5_ERROR_IF_USED` (`ghl_error`), terminating the process.
- **Test coverage:** `unit_test_grhayl_core_test_suite.c` tests simple/hybrid bounds, rectangular extrema, and constant atmosphere; its tabulated setup/atmosphere checks remain commented out (`TODO`). Error paths live in `unit_test_code_error.c`.

Evidence (descriptive behavior): `GRHayL/GRHayL_Core/initialize_eos.c::ghl_initialize_tabulated_eos`; corroboration: `GRHayL/include/ghl_nrpyeos_tabulated.h::GHL_HDF5_ERROR_IF_USED`.

## Sources

- [GRHayL/GRHayL_Core/initialize_eos.c](../../GRHayL/GRHayL_Core/initialize_eos.c) - EOS wrappers, global pointer storage, and dispatch
- [GRHayL/include/ghl.h](../../GRHayL/include/ghl.h) - `ghl_eos_parameters`, `ghl_error_used_disabled_hdf5`
- [GRHayL/include/ghl_eos_functions.h](../../GRHayL/include/ghl_eos_functions.h) - `extern` EOS function-pointer declarations
- [GRHayL/include/ghl_eos_functions_declaration.h](../../GRHayL/include/ghl_eos_functions_declaration.h) - internal storage definitions
- [GRHayL/include/ghl_nrpyeos_tabulated.h](../../GRHayL/include/ghl_nrpyeos_tabulated.h) - `GHL_HDF5_ERROR_IF_USED`
- [GRHayL/EOS/Hybrid/NRPyEOS_initialize_hybrid_functions.c](../../GRHayL/EOS/Hybrid/NRPyEOS_initialize_hybrid_functions.c) - hybrid pointer assignment
- [GRHayL/EOS/Tabulated/NRPyEOS_initialize_tabulated_functions.c](../../GRHayL/EOS/Tabulated/NRPyEOS_initialize_tabulated_functions.c) - tabulated pointer assignment
- [Unit_Tests/unit_test_grhayl_core_test_suite.c](../../Unit_Tests/unit_test_grhayl_core_test_suite.c) - simple/hybrid bounds and commented tabulated checks
- [Unit_Tests/unit_test_code_error.c](../../Unit_Tests/unit_test_code_error.c) - invalid EOS parameter error keys

## See Also

- Implements: [EOS Initialization And Dispatch](../gems/eos/initialization-and-dispatch.md)
- Depends on: [Build And Configure](../architecture/build-and-configure.md)
- See also: [Atmosphere](../gems/atmosphere.md)
- Validated by: [Core Tests And Reference Data](tests-and-reference-data.md)
