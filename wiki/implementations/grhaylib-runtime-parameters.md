# GRHayLib Runtime Parameters

> Global allocation/cleanup, `GRHayLib_paramcheck`, Cactus-to-GRHayL parameter mapping, EOS branches, Con2Prim keyword parsing, PPM overrides, and steering limits in `implementations/GRHayLib/`. · Status: confirmed
> Up: [Implementations](index.md)

## Summary

`implementations/GRHayLib/src/initialize_and_shutdown.c` allocates and releases the thorn's shared global pointers `ghl_params` and `ghl_eos`, validates Cactus parameters from `param.ccl`, translates Con2Prim and EOS keywords, overrides Core PPM defaults, and initializes simple, hybrid, or tabulated EOS dispatch. Cactus `STEERABLE=RECOVER` parameters do not auto-refresh `ghl_params` or `ghl_eos`.

## Detail

### Global Allocation, Initialization, And Termination

`GRHayLib_initialize` allocates zeroed `ghl_params` and `ghl_eos` via checked `calloc` (freeing partial storage and clearing globals on allocation failure), parses `con2prim_routine` and `con2prim_backup_routines[3]` via `parse_C2P_routine_keyword`, calls `GRHayLib_paramcheck` and `ghl_initialize_params`, copies PPM Cactus parameters into `ghl_params`, and initializes EOS parameters and dispatch. Later fatal setup errors have no GRHayLib-local rollback; tabulated initialization still leaves a tagged empty aggregate on failure ([Tabulated Table Layout](../gems/eos/tabulated-table-layout.md)).

`GRHayLib_terminate` checks `ghl_eos`: if non-null and `ghl_eos->eos_type == ghl_eos_tabulated`, it calls `ghl_tabulated_free_memory(ghl_eos)`, frees `ghl_eos` and `ghl_params`, and sets both to `NULL`. Null globals and repeated calls after cleanup are safe; non-null globals must denote valid heap allocations with initialized cleanup dispatch. When `ID_converter_ILGRMHD` is active, `schedule.ccl` skips `GRHayLib_initialize` while still scheduling `GRHayLib_terminate`, which releases globals allocated by GRHayLET's downstream `IllinoisGRMHD_backward_compatible_initialize` (which has no termination hook of its own).

Evidence (descriptive behavior): `initialize_and_shutdown.c::GRHayLib_initialize`; corroboration: `initialize_and_shutdown.c::GRHayLib_terminate` and `schedule.ccl`.

### Parameter Mapping, PPM Overrides, And Steering Limits

`GRHayLib_initialize` passes parsed main/backup Con2Prim IDs, `evolve_entropy`, `evolve_temperature`, `calc_primitive_guess`, `Psi6threshold`, `max_Lorentz_factor`, and `Lorenz_damping_factor` into `ghl_initialize_params` ([Shared Parameters And Enums](../core/shared-parameters-and-enums.md)). `GRHayLib_paramcheck` requires finite `max_Lorentz_factor >= 1` (for subluminal velocity-rescaling edge cases and retry bounds at extreme Lorentz factors, see [Velocity And `u0`](../core/velocity-u0.md)). After `ghl_initialize_params`, GRHayLib overwrites Core PPM defaults with `ppm_flattening_epsilon`, `ppm_flattening_omega1`, `ppm_flattening_omega2`, `ppm_shock_epsilon`, `ppm_shock_eta1`, `ppm_shock_eta2`, and `ppm_shock_k0` ([PPM Flow](../gems/reconstruction/ppm-flow.md)).

Although `param.ccl` marks selected solver, bound, table-path, and PPM parameters `STEERABLE=RECOVER`, `schedule.ccl` provides no post-recovery refresh hook and re-invoking `GRHayLib_initialize` would overwrite live pointers without freeing them (`doc/documentation.tex`). Any downstream runtime update must maintain derived fields, free tabulated/NN allocations before reinitialization, and keep global EOS function pointers synchronized.

Evidence (descriptive behavior): `initialize_and_shutdown.c::GRHayLib_paramcheck`; corroboration: `param.ccl::STEERABLE=RECOVER`.

### EOS Branches And Required/Optional Limits

| `EOS_type` | Required inputs & optional `-1` sentinels | Dispatch & initialization calls |
| --- | --- | --- |
| `"Simple"` | Finite `Gamma > 1`, `rho_b_atm`, `P_atm`; optional `rho_b_{min,max}`, `P_{min,max}` (Core defaults `0.0` / `1e300`). (Ev: `initialize_and_shutdown.c::ghl_initialize_simple_eos_functions_and_params`) | Sets `ghl_con2prim_multi_method = ghl_con2prim_hybrid_multi_method`; calls `ghl_initialize_simple_eos_functions_and_params` (installs `ghl_eos_simple` dispatch). |
| `"Hybrid"` | `Gamma_th`, `Gamma_ppoly_in[0]`, `k_ppoly0`, `rho_b_atm`, plus `neos - 1` breakpoints and `neos` gammas; optional `rho_b_{min,max}` (disabled `P`/`eps`/`S` bounds map to `+/-DBL_MAX`). (Ev: `initialize_and_shutdown.c::ghl_initialize_hybrid_eos_functions_and_params`) | Sets `ghl_con2prim_multi_method = ghl_con2prim_hybrid_multi_method`; calls `ghl_initialize_hybrid_eos_functions_and_params` (installs `ghl_eos_hybrid` dispatch). |
| `"Tabulated"` | `rho_b_atm`, `Y_e_atm`, `T_atm`, nonempty `EOS_tablepath`; optional `rho_b`/`Y_e`/`T` bounds clamped to table limits. (Ev: `initialize_and_shutdown.c::ghl_initialize_tabulated_eos`) | Sets `ghl_con2prim_tabulated_multi_method`, calls `ghl_initialize_eos_functions(ghl_eos_tabulated)`, then `ghl_initialize_tabulated_eos` with `parse_eos_table_type_keyword` (`"stellarcollapse"` -> `ghl_eos_table_stellarcollapse`, else `ghl_eos_table_unknown`), `eos_table_clean_sound_speed`, and `enable_backup_nn_primitive_guess` as `enable_neural_net_c2p`. |

Empty or unrecognized `EOS_type` errors. Core still rejects several inputs that `param.ccl` / `GRHayLib_paramcheck` admit (`rho_b_atm = 0` for simple/hybrid, used `Gamma_ppoly_in` in `{0, 1}`, `Gamma_th = 1`, nonpositive/nonincreasing `rho_ppoly_in`), halting via `ghl_abort_if_error`. Simple/hybrid Core initialization also sets unused `Y_e_atm` and `T_atm` placeholders to zero before constant-atmosphere copies.

### Con2Prim Keyword Parsing And Tabulated NN Guess

`parse_C2P_routine_keyword` maps `Noble2D`, `Noble1D`, `Noble1D_entropy`, `Noble1D_entropy2`, `Font1D`, `Palenzuela1D`, `Palenzuela1D_entropy`, `Newman1D`, `Newman1D_entropy`, and `None` (exposed in `param.ccl` only for `con2prim_backup_routines[3]`) to `ghl_con2prim_id_t`, returning `-100` for unknown strings. `GRHayLib_paramcheck`:

- rejects EOS-incompatible solver choices per EOS family (Simple rejects `Newman1D`/`Font1D`/`Newman1D_entropy`; Hybrid rejects `Newman1D`/`Newman1D_entropy`; Tabulated rejects `Noble1D`/`Noble1D_entropy`/`Noble1D_entropy2`/`Font1D`; see [Solver Matrix](../gems/con2prim/solver-matrix.md));
- reports rejected backup choices with their exact slot number and value;
- rejects all entropy methods (`Noble1D_entropy`, `Noble1D_entropy2`, `Palenzuela1D_entropy`, `Newman1D_entropy`) when `evolve_entropy` is false;
- restricts `Noble1D_entropy` and `Noble1D_entropy2` to Simple EOS or equivalent one-piece Hybrid metadata (`Gamma_th == Gamma_ppoly_in[0]`).

Because ordinary hybrid Noble recovery can reject marginal cold-boundary states after velocity-limit finalization exposes nonpositive pressure, deployments needing those states should configure a robust backup such as `Font1D`. `enable_backup_nn_primitive_guess` is not a `ghl_con2prim_id_t` solver: it sets `enable_neural_net_c2p` on tabulated EOS initialization so `ghl_con2prim_tabulated_multi_method` can retry via `ghl_c2p_nn_guess_primitives` after main/backup failure ([Neural-Network Primitive Guess](../gems/con2prim/neural-network-primitive-guess.md)).

Evidence (descriptive behavior): `initialize_and_shutdown.c::parse_C2P_routine_keyword`; corroboration: `param.ccl::enable_backup_nn_primitive_guess`.

## Sources

- [implementations/GRHayLib/src/initialize_and_shutdown.c](../../implementations/GRHayLib/src/initialize_and_shutdown.c) - `GRHayLib_paramcheck`, `GRHayLib_initialize`, `GRHayLib_terminate`, `parse_C2P_routine_keyword`, `parse_eos_table_type_keyword`
- [implementations/GRHayLib/param.ccl](../../implementations/GRHayLib/param.ccl) - Cactus parameter definitions, ranges, and `STEERABLE=RECOVER`
- [implementations/GRHayLib/schedule.ccl](../../implementations/GRHayLib/schedule.ccl) - `ID_converter_ILGRMHD` guard and `CCTK_TERMINATE`
- [implementations/GRHayLib/doc/documentation.tex](../../implementations/GRHayLib/doc/documentation.tex) - steering caveats
- [GRHayLET backward_compatible_initialize.c](https://github.com/GRHayL/GRHayLET/blob/main/IllinoisGRMHD/src/backward_compatible_initialize.c) - downstream compatibility allocation and Hybrid EOS setup
- [GRHayLET IllinoisGRMHD schedule.ccl](https://github.com/GRHayL/GRHayLET/blob/main/IllinoisGRMHD/schedule.ccl) - downstream schedule under `ID_converter_ILGRMHD`

## See Also

- Depends on: [GRHayLib Cactus Build](grhaylib-cactus-build.md)
- Depends on: [Solver Matrix](../gems/con2prim/solver-matrix.md)
- Depends on: [EOS Initialization And Dispatch](../gems/eos/initialization-and-dispatch.md)
- See also: [Neural-Network Primitive Guess](../gems/con2prim/neural-network-primitive-guess.md)
- See also: [GRHayLib Verification](grhaylib-verification.md)
