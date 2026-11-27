# Con2Prim Recovery Flow

> Seven-step runtime recovery sequence, densitized versus undensitized boundaries, backup/NN retry rules, and `GHL_DISABLE_HDF5` stubs. · Status: confirmed
> Up: [Con2Prim](index.md)

## Summary

A logical Con2Prim recovery proceeds from diagnostics initialization and optional conservative limiting on densitized variables, through undensitization, optional primitive guessing, main/backup (and optional tabulated NN) solver retries on undensitized conservatives, to post-recovery primitive limiting/`u0` and conservative/stress-energy recomputation. Failed recoveries may leave partial primitive mutations and do not automatically reset to atmosphere.

## Detail

### Ordered Recovery Sequence And Densitized Boundary

1. **Initialize diagnostics:** `ghl_initialize_diagnostics` (`initialize_diagnostics.c`) zeroes `n_iter`, clears `tau_fix`, `Stilde_fix`, `speed_limited`, `nn_guess_used`, and `backup[0..2]`, and sets `which_routine = ghl_con2prim_id_None`.
2. **Apply conservative limits (densitized):** `ghl_apply_conservative_limits` (`apply_conservative_limits.c`) mutates densitized `tau` and `SD` using `prims->BU` and updates `tau_fix` / `Stilde_fix`.
3. **Undensitize conservatives:** `ghl_undensitize_conservatives(metric_adm.sqrt_detgamma, &cons, &cons_undens)` (`undensitize_conservatives.c`); all selectors, drivers, and `ghl_guess_primitives` expect `cons_undens`.
4. **Optional primitive guess:** when `params->calc_prim_guess` is true, multi-method drivers invoke `ghl_guess_primitives` before the main solver (hybrid/simple cold pressure; tabulated Palenzuela-style estimate). When skipped, callers must supply solver-required inputs including `Y_e` and `temperature` for tabulated Palenzuela.
5. **Main solve & up to three backups:** single-shot selectors (`ghl_con2prim_{hybrid,tabulated}_select_method`) run one `c2p_key`. Multi-method drivers save the initial guess, run `params->main_routine`, and on failure iterate `params->backup_routine[0..2]`, stopping at success, slot exhaustion, or `ghl_con2prim_id_None`. Before each backup, the driver sets `diagnostics->backup[n] = true` and restores `*prims` to the saved guess. Attempt flags (`speed_limited`, `backup[]`, `nn_guess_used`) stay sticky across failed attempts; on final failure `which_routine` stays `None`, `n_iter` is unspecified, and `*prims` is not rolled back.
6. **Post-recovery limits & `u0`:** `ghl_enforce_primitive_limits_and_compute_u0` clamps EOS-specific thermodynamic bounds and runs `ghl_limit_v_and_compute_u0` only if thermodynamic updates succeed.
7. **Recompute densitized conservatives (and `Tmunu`):** `ghl_compute_conservs` or `ghl_compute_conservs_and_Tmunu` consumes valid `prims->eps` and `prims->u0` to produce consistent densitized `cons` (and `Tmunu`).

| Routine | Conservative input | Conservative output | Boundary note |
| --- | --- | --- | --- |
| `ghl_apply_conservative_limits` | Densitized `cons` (`prims->BU` only) | Mutated densitized `cons` | Updates `tau_fix` and `Stilde_fix`. (Ev: `GRHayL/Con2Prim/apply_conservative_limits.c::ghl_apply_conservative_limits`) |
| `ghl_undensitize_conservatives` | Densitized `cons`, `psi6` | Undensitized `cons_undens` | Divides `rho`, `SD`, `tau`, `Y_e`, `entropy` by `psi6`. (Ev: `GRHayL/Con2Prim/undensitize_conservatives.c::ghl_undensitize_conservatives`) |
| `ghl_guess_primitives` & selectors/drivers | Undensitized `cons_undens` | Primitives only | Multi-method saves/restores guess across backups. (Ev: `GRHayL/Con2Prim/con2prim_multi_method.c::ghl_con2prim_hybrid_multi_method`) |
| `ghl_enforce_primitive_limits_and_compute_u0` | None | None | Recomputes EOS fields and `u0`; failure can leave partial writes. (Ev: `GRHayL/Con2Prim/enforce_primitive_limits_and_compute_u0.c`) |
| `ghl_compute_conservs` / `_and_Tmunu` | Primitives (`eps`, `u0` required) | Densitized `cons` (+ `Tmunu`) | Post-limit synchronization. (Ev: `GRHayL/Con2Prim/compute_conservs.c::ghl_compute_conservs`) |

### Tabulated NN Retry Hook And `GHL_DISABLE_HDF5` Stubs

- **Tabulated NN retry:** Separate from `params->calc_prim_guess`, when `eos->enable_neural_net_c2p` is true and the first tabulated main attempt fails, `ghl_con2prim_tabulated_multi_method` builds `prims_guess_nn` via `ghl_c2p_nn_guess_primitives`, sets `diagnostics->nn_guess_used = true`, and retries the main routine; if later backup slots fail with the ordinary guess, the driver reuses `prims_guess_nn` to retry each failed backup. `nn_guess_used` records an attempted NN retry, not success ([Neural-Network Primitive Guess](neural-network-primitive-guess.md)).
- **`GHL_DISABLE_HDF5`:** Tabulated select/multi-method, the tabulated branch of `ghl_enforce_primitive_limits_and_compute_u0`, and the link-visible non-mutating direct tabulated solver stubs in `GRHayL/Con2Prim/hdf5_disabled_stubs.c` all return `ghl_error_used_disabled_hdf5` before dispatching.

Evidence (descriptive behavior): `GRHayL/Con2Prim/con2prim_multi_method.c::ghl_con2prim_tabulated_multi_method`; corroboration: `GRHayL/Con2Prim/hdf5_disabled_stubs.c`.

## Sources

- [GRHayL/include/ghl_con2prim.h](../../../GRHayL/include/ghl_con2prim.h) - recovery, diagnostics, and helper declarations
- [GRHayL/Con2Prim/con2prim_multi_method.c](../../../GRHayL/Con2Prim/con2prim_multi_method.c) - selectors, backup loop, NN retry hook, and no-HDF5 guards
- [GRHayL/Con2Prim/initialize_diagnostics.c](../../../GRHayL/Con2Prim/initialize_diagnostics.c) - diagnostics reset
- [GRHayL/Con2Prim/apply_conservative_limits.c](../../../GRHayL/Con2Prim/apply_conservative_limits.c) - pre-recovery conservative limits
- [GRHayL/Con2Prim/undensitize_conservatives.c](../../../GRHayL/Con2Prim/undensitize_conservatives.c) - undensitization helper
- [GRHayL/Con2Prim/guess_primitives.c](../../../GRHayL/Con2Prim/guess_primitives.c) - ordinary primitive guess
- [GRHayL/Con2Prim/enforce_primitive_limits_and_compute_u0.c](../../../GRHayL/Con2Prim/enforce_primitive_limits_and_compute_u0.c) - post-recovery limits and `u0`
- [GRHayL/Con2Prim/compute_conservs.c](../../../GRHayL/Con2Prim/compute_conservs.c) - densitized conservative recompute
- [GRHayL/Con2Prim/hdf5_disabled_stubs.c](../../../GRHayL/Con2Prim/hdf5_disabled_stubs.c) - non-mutating no-HDF5 direct tabulated solver stubs

## See Also

- Depends on: [Solver Matrix](solver-matrix.md)
- Depends on: [Limits And Conversions](limits-and-conversions.md)
- See also: [Neural-Network Primitive Guess](neural-network-primitive-guess.md)
- Validated by: [Con2Prim Tests And Reference Data](tests-and-reference-data.md)
