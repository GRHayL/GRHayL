# EOS Initialization And Dispatch

> Simple, hybrid, and tabulated `ghl_eos_parameters` setup, caller validation boundaries, dispatch outcomes, and `GHL_DISABLE_HDF5` behavior. · Status: confirmed
> Up: [Equation Of State (EOS)](index.md)

## Summary

GRHayL initializes `ghl_eos_parameters` (`ghl_eos_simple`, `ghl_eos_hybrid`, `ghl_eos_tabulated`) via public `*_functions_and_params` wrappers or low-level parameter initializers backed by process-wide function-pointer globals ([Core EOS Dispatch](../../core/eos-dispatch.md)). Simple and hybrid setup construct local candidates leaving destinations untouched on failure, whereas tabulated setup transfers heap allocations and requires freeing live tables before reinitialization.

## Detail

### Family Parameter Setup And Validation Rules

- **Simple EOS (`ghl_initialize_simple_eos`):** uses one-piece ideal-fluid setup (`neos = 1`, `Gamma_th = Gamma_ppoly[0] = Gamma`), normalizes negative density/pressure min to `0` and max to `1e300`, derives rectangular `eps`/`entropy` extrema from opposite density endpoints for `Gamma > 1` (`DBL_MAX` at zero density floor or `1e300` pressure ceiling; `0` at zero pressure ceiling), sets `tau_atm = rho_atm * eps_atm`, and zeroes unused `Y_e_atm` and `T_atm`.
- **Hybrid EOS (`ghl_initialize_hybrid_eos`):** requires `1 <= neos <= MAX_EOS_PARAMS`, finite nonsingular gammas (`Gamma_ppoly` not in `{0, 1}`, `Gamma_th != 1`), and `neos - 1` finite positive strictly increasing density breakpoints (`NULL` permitted for one piece). Computes `K_ppoly`, `eps_integ_const`, and `neos - 1` transitions in `p_ppoly` (leaving `p_ppoly[neos - 1] = 0`); zero `K_ppoly0` intentionally gives a zero cold-pressure curve (with no unique cold density inverse); otherwise every `K_ppoly` must be nonzero. Uses `-DBL_MAX`/`DBL_MAX` for unbounded pressure/energy/entropy extrema and zeroes `Y_e_atm`, `T_atm`.
- **Tabulated EOS (`ghl_initialize_tabulated_eos`):** records `table_type` and `clean_sound_speed`, nulls beta-equilibrium and NN pointers before fallible I/O, rejects nonfinite inputs before reading, clamps requested min/max bounds to table limits (**atmosphere `(rho, Y_e, T)` is not clamped** and must lie inside table bounds), computes atmosphere `press`/`eps`/`entropy`, sets `tau_atm = rho_min * eps_min`, sets `root_finding_precision = 1e-10`, and checks finiteness of all published bounds and extrema. Public `ghl_initialize_tabulated_eos_functions_and_params` hard-codes `ghl_eos_table_stellarcollapse`, `clean_sound_speed = false`, and `enable_neural_net_c2p = false`.

Evidence (descriptive behavior): `GRHayL/GRHayL_Core/initialize_eos.c::ghl_initialize_hybrid_eos`; corroboration: `GRHayL/GRHayL_Core/initialize_eos.c::ghl_initialize_simple_eos` and `GRHayL/GRHayL_Core/initialize_eos.c::ghl_initialize_tabulated_eos`.

### Dispatch Outcomes And `GHL_DISABLE_HDF5`

For process-wide pointer storage (`ghl_eos_functions_declaration.h`), stage ordering, and no-rollback behavior, see [Core EOS Dispatch](../../core/eos-dispatch.md). In summary:
- `NRPyEOS_initialize_hybrid_functions` wires hybrid pointers and `ghl_compute_h_and_cs2`; `NRPyEOS_initialize_tabulated_functions` (HDF5 builds) wires tabulated pointers, beta-equilibrium cleanup (`ghl_tabulated_free_beq_quantities`), and tabulated `ghl_compute_h_and_cs2`.
- `ghl_initialize_eos_functions` routes `ghl_compute_h_and_cs2` and `ghl_con2prim_multi_method` to hybrid (simple/hybrid) or tabulated (HDF5 builds); unknown `ghl_eos_t` returns `void` with partial dispatch state. Deprecated generic `ghl_calculate_HLLE_fluxes_dirn*` globals are never assigned.
- Under `GHL_DISABLE_HDF5`, tabulated parameter initializers return `ghl_error_used_disabled_hdf5` (publishing an empty `ghl_eos_tabulated` aggregate on valid storage), while `ghl_initialize_eos_functions(ghl_eos_tabulated)` terminates via `GHL_HDF5_ERROR_IF_USED`.

Evidence (descriptive behavior): `GRHayL/EOS/Hybrid/NRPyEOS_initialize_hybrid_functions.c`; corroboration: `GRHayL/EOS/Tabulated/NRPyEOS_initialize_tabulated_functions.c`.

## Sources

- [GRHayL/GRHayL_Core/initialize_eos.c](../../../GRHayL/GRHayL_Core/initialize_eos.c) - simple, hybrid, and tabulated initializers and dispatch selection
- [GRHayL/include/ghl.h](../../../GRHayL/include/ghl.h) - `ghl_eos_parameters`, `ghl_eos_t`, `ghl_eos_table_t`
- [GRHayL/include/ghl_eos_functions.h](../../../GRHayL/include/ghl_eos_functions.h) - public EOS function pointers
- [GRHayL/EOS/Hybrid/NRPyEOS_initialize_hybrid_functions.c](../../../GRHayL/EOS/Hybrid/NRPyEOS_initialize_hybrid_functions.c) - hybrid function-pointer setup
- [GRHayL/EOS/Tabulated/NRPyEOS_initialize_tabulated_functions.c](../../../GRHayL/EOS/Tabulated/NRPyEOS_initialize_tabulated_functions.c) - tabulated function-pointer setup
- [docs/raw/EOS.dox](../../../docs/raw/EOS.dox) - EOS overview documentation

## See Also

- Depends on: [EOS Dispatch](../../core/eos-dispatch.md)
- Implements: [Hybrid Piecewise-Polytrope EOS](hybrid-piecewise-polytrope.md)
- Implements: [Tabulated Table Layout](tabulated-table-layout.md)
- See also: [Tabulated Interpolator Catalog](tabulated-interpolator-catalog.md)
