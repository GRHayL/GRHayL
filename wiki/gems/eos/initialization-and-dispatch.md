# EOS Initialization And Dispatch

> Simple, hybrid, and tabulated `ghl_eos_parameters` setup, caller validation rules, rectangular extrema, and `tau_atm`. · Status: confirmed
> Up: [Equation Of State (EOS)](index.md)

## Summary

GRHayL initializes `ghl_eos_parameters` (`ghl_eos_simple`, `ghl_eos_hybrid`, `ghl_eos_tabulated`) via public `*_functions_and_params` wrappers or low-level parameter initializers backed by process-wide function-pointer globals ([Core EOS Dispatch](../../core/eos-dispatch.md)). Simple and hybrid setup construct local candidates leaving destinations untouched on failure, whereas tabulated setup transfers heap allocations and requires freeing live tables before reinitialization.

## Detail

### Family Parameter Setup And Validation Rules

- **Simple EOS (`ghl_initialize_simple_eos`):** rejects `NULL` output, nonfinite scalars, `rho_atm <= 0`, `press_atm < 0`, `Gamma` in `{0, 1}`, `rho_max == 0`, and inverted `rho_min > rho_max` or `press_min > press_max`; uses one-piece ideal-fluid setup (`neos = 1`, `Gamma_th = Gamma_ppoly[0] = Gamma`), normalizes negative density/pressure min to `0` and max to `1e300`, derives rectangular `eps`/`entropy` extrema by pairing `press_min` with `rho_max` and `press_max` with `rho_min` for `Gamma > 1` (`DBL_MAX` at zero `rho_min` or `1e300` ceiling; `0` at zero ceiling), sets `tau_atm = rho_atm * eps_atm`, and zeroes unused `Y_e_atm` and `T_atm`.
- **Hybrid EOS (`ghl_initialize_hybrid_eos`):** rejects `NULL` output, nonfinite scalars, `rho_atm <= 0`, `rho_max == 0`, and `rho_min > rho_max`, and requires `1 <= neos <= MAX_EOS_PARAMS`, finite nonsingular gammas (`Gamma_ppoly` not in `{0, 1}`, `Gamma_th != 1`), and `neos - 1` finite positive strictly increasing density breakpoints (`NULL` permitted for one piece). Computes `K_ppoly`, `eps_integ_const`, and `neos - 1` transitions in `p_ppoly` (leaving `p_ppoly[neos - 1] = 0`); zero `K_ppoly0` intentionally gives a zero cold-pressure curve (with no unique cold density inverse); otherwise every `K_ppoly` must be nonzero. Uses `-DBL_MAX`/`DBL_MAX` for unbounded pressure/energy/entropy extrema, sets `tau_atm = rho_atm * eps_atm`, and zeroes `Y_e_atm`, `T_atm`.
- **Tabulated EOS (`ghl_initialize_tabulated_eos`):** records `table_type` and `clean_sound_speed`, nulls beta-equilibrium and NN pointers before fallible I/O, rejects nonfinite inputs before reading, clamps requested min/max bounds to table limits (**atmosphere `(rho, Y_e, T)` is not clamped** and must lie inside table bounds), computes atmosphere `press`/`eps`/`entropy`, sets `tau_atm = rho_min * eps_min`, sets `root_finding_precision = 1e-10`, and checks finiteness of all published bounds and extrema. Public `ghl_initialize_tabulated_eos_functions_and_params` hard-codes `ghl_eos_table_stellarcollapse`, `clean_sound_speed = false`, and `enable_neural_net_c2p = false`.

Evidence (descriptive behavior): `GRHayL/GRHayL_Core/initialize_eos.c::ghl_initialize_hybrid_eos`; corroboration: `GRHayL/GRHayL_Core/initialize_eos.c::ghl_initialize_simple_eos` and `GRHayL/GRHayL_Core/initialize_eos.c::ghl_initialize_tabulated_eos`.

For process-wide function-pointer storage (`ghl_eos_functions_declaration.h`), `ghl_initialize_eos_functions` stage ordering, no-rollback wrapper behavior, deprecated `ghl_calculate_HLLE_fluxes_dirn*` globals, and `GHL_DISABLE_HDF5` dispatch divergence, see [Core EOS Dispatch](../../core/eos-dispatch.md).

## Sources

- [GRHayL/GRHayL_Core/initialize_eos.c](../../../GRHayL/GRHayL_Core/initialize_eos.c) - simple, hybrid, and tabulated initializers and validation rules
- [GRHayL/include/ghl.h](../../../GRHayL/include/ghl.h) - `ghl_eos_parameters`, `ghl_eos_t`, `ghl_eos_table_t`
- [docs/raw/EOS.dox](../../../docs/raw/EOS.dox) - EOS overview documentation

## See Also

- Depends on: [EOS Dispatch](../../core/eos-dispatch.md)
- Implements: [Hybrid Piecewise-Polytrope EOS](hybrid-piecewise-polytrope.md)
- Implements: [Tabulated Table Layout](tabulated-table-layout.md)
- See also: [Tabulated Interpolator Catalog](tabulated-interpolator-catalog.md)
