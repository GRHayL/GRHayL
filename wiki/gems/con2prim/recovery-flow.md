# Con2Prim Recovery Flow

This page summarizes runtime recovery order and conservative-variable data
boundaries for `GRHayL/Con2Prim/`. It is a route into source, Doxygen, and
tests; repo-local source remains authority.

## Ordered Flow

1. Initialize `ghl_con2prim_diagnostics` with
   `ghl_initialize_diagnostics`. It clears `tau_fix`, `Stilde_fix`,
   `speed_limited`, all three `backup` slots, and sets `which_routine` to
   `ghl_con2prim_id_None`; the diagnostics struct and public helper list are
   declared in `GRHayL/include/ghl_con2prim.h`.
   Source: `GRHayL/Con2Prim/initialize_diagnostics.c`.
2. Apply conservative limits when the caller has densitized conservative data
   that needs the Faber-style inequalities before recovery. The helper mutates
   densitized `tau` and `SD`, needs only `prims->BU` from primitives, and
   updates `diagnostics->tau_fix` or `diagnostics->Stilde_fix`.
   Source: `GRHayL/Con2Prim/apply_conservative_limits.c`;
   test: `Unit_Tests/unit_test_apply_conservative_limits.c`.
3. Convert densitized conservative fields to undensitized fields with
   `ghl_undensitize_conservatives(metric_adm.sqrt_detgamma, &cons,
   &cons_undens)`. The Con2Prim drivers and solvers take the undensitized
   object, not evolved densitized storage.
   Source: `GRHayL/Con2Prim/undensitize_conservatives.c`;
   tests: `Unit_Tests/unit_test_con2prim_multi_method_hybrid.c`,
   `Unit_Tests/unit_test_con2prim_tabulated.c`.
4. Optionally guess primitives. `params->calc_prim_guess` makes both
   multi-method drivers call `ghl_guess_primitives` before the main attempt.
   The helper contract takes undensitized conservatives; hybrid/simple guesses
   set cold-pressure data, while tabulated guesses use a Palenzuela-style
   estimate and table bounds.
   Source: `GRHayL/Con2Prim/guess_primitives.c`,
   `GRHayL/Con2Prim/con2prim_multi_method.c`.
5. Call either a selector directly or the multi-method driver. Selectors
   (`ghl_con2prim_hybrid_select_method`,
   `ghl_con2prim_tabulated_select_method`) dispatch one `c2p_key` and return
   its error. Multi-method drivers choose `params->main_routine`, own optional
   guessing, and own backup retries.
   Source: `GRHayL/Con2Prim/con2prim_multi_method.c`.
6. If the main routine fails inside a multi-method driver, try up to three
   `params->backup_routine` slots. The loop stops after success, after three
   slots, or immediately at the first `ghl_con2prim_id_None`. For each backup
   actually tried, it sets `diagnostics->backup[n] = true`, resets `*prims` to
   the saved primitive guess, and calls the same selector layer again. The
   successful solver updates `diagnostics->which_routine`; solver bodies also
   update `diagnostics->n_iter` where applicable.
   Source: `GRHayL/Con2Prim/con2prim_multi_method.c`;
   tests: `Unit_Tests/unit_test_hybrid_failure.c`,
   `Unit_Tests/unit_test_con2prim_tabulated.c`.
7. Enforce primitive limits and compute `u0` with
   `ghl_enforce_primitive_limits_and_compute_u0`. The helper applies EOS-
   specific primitive floors/ceilings, recomputes thermodynamic fields, then
   applies the speed limiter and computes `u0`.
   Source: `GRHayL/Con2Prim/enforce_primitive_limits_and_compute_u0.c`;
   tests: `Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c`,
   `Unit_Tests/data_gen/unit_test_data_con2prim_multi_method_hybrid.c`.
8. Recompute conservatives, and stress-energy when the caller needs it, after
   primitive limiting. `ghl_compute_conservs` returns densitized conservatives;
   `ghl_compute_conservs_and_Tmunu` returns densitized conservatives plus
   `Tmunu`. Both require initialized `prims->eps` and `prims->u0`.
   Source: `GRHayL/Con2Prim/compute_conservs.c`,
   `GRHayL/Con2Prim/compute_conservs_and_Tmunu.c`;
   tests: `Unit_Tests/unit_test_compute_conservs_and_Tmunu.c`,
   `Unit_Tests/data_gen/unit_test_data_con2prim_multi_method_hybrid.c`.

## Tabulated NN Retry Hook

Tabulated neural-network primitive guesses route through
[neural-network primitive guess](neural-network-primitive-guess.md). They do
not replace `params->calc_prim_guess`: that flag controls the ordinary
`ghl_guess_primitives` call before the first main routine attempt.
`eos->enable_neural_net_c2p` instead controls fallback retries for tabulated
recovery. After the first failed tabulated main attempt, the driver computes
`prims_guess_nn` with `ghl_c2p_nn_guess_primitives`, resets `*prims` to that
guess, and retries the main selector. If later tabulated backup attempts fail,
the driver reuses the stored `prims_guess_nn` and calls the same backup
selector again.

Sources: `GRHayL/Con2Prim/con2prim_multi_method.c`,
`GRHayL/Con2Prim/Tabulated/neural_network_guess/c2p_nn_guess_primitives.c`,
`GRHayL/include/ghl.h`, `GRHayL/include/ghl_con2prim.h`.
Tests: `Unit_Tests/unit_test_c2p_nn_guess.c`,
`Unit_Tests/unit_test_con2prim_tabulated.c`.

## Selector And Multi-Method Split

- Selectors are single-shot dispatchers. They map one `ghl_con2prim_id_t` to
  an EOS-compatible solver or return `ghl_error_invalid_c2p_key`.
  Source: `GRHayL/Con2Prim/con2prim_multi_method.c`.
- Multi-method drivers are policy wrappers. They read `params->main_routine`,
  optionally compute a primitive guess, preserve that guess, and run the backup
  loop over the three configured slots.
  Source: `GRHayL/include/ghl.h`,
  `GRHayL/GRHayL_Core/initialize_params.c`,
  `GRHayL/Con2Prim/con2prim_multi_method.c`.
- Diagnostics lifecycle is split. Callers initialize diagnostics; conservative
  limits record `tau_fix` and `Stilde_fix`; multi-method records backup slots;
  successful solver implementations set `which_routine` and often `n_iter`;
  post-limits write `speed_limited`.
  Source: `GRHayL/include/ghl_con2prim.h`,
  `GRHayL/Con2Prim/initialize_diagnostics.c`,
  `GRHayL/Con2Prim/apply_conservative_limits.c`,
  `GRHayL/Con2Prim/enforce_primitive_limits_and_compute_u0.c`,
  solver files under `GRHayL/Con2Prim/Hybrid/` and
  `GRHayL/Con2Prim/Tabulated/`.

## HDF5-Disabled Tabulated Behavior

When built with `GHL_DISABLE_HDF5`, tabulated select and tabulated
multi-method return `ghl_error_used_disabled_hdf5` before solver dispatch or
backup behavior. The tabulated branch in
`ghl_enforce_primitive_limits_and_compute_u0` returns the same error. The
code-error tests skip HDF5-only cases in no-HDF5 builds.

Sources: `GRHayL/Con2Prim/con2prim_multi_method.c`,
`GRHayL/Con2Prim/enforce_primitive_limits_and_compute_u0.c`,
`GRHayL/include/ghl.h`, `Unit_Tests/unit_test_code_error.c`.

## Densitized Boundary Table

| Function or driver | Conservative input | Conservative output | Boundary contract |
| --- | --- | --- | --- |
| `ghl_apply_conservative_limits` | Densitized `cons`; only `prims->BU` required from primitives | Mutates densitized `cons` | Pre-recovery limiter for evolved `tau` and `SD`; records `tau_fix` and `Stilde_fix`. Source: `GRHayL/Con2Prim/apply_conservative_limits.c`; test: `Unit_Tests/unit_test_apply_conservative_limits.c`. |
| `ghl_undensitize_conservatives` | Densitized `cons` plus `psi6` | Undensitized `cons_undens` | Divides `rho`, `SD`, `tau`, `Y_e`, and `entropy` by `psi6`; use before Con2Prim solvers. Source: `GRHayL/Con2Prim/undensitize_conservatives.c`. |
| `ghl_guess_primitives` | Undensitized `cons_undens` by source contract | Primitive guess only | Multi-method drivers pass undensitized data. Some direct hybrid fixture code calls the helper separately while also calling the selector; use the source/helper contract for new calls. Source: `GRHayL/Con2Prim/guess_primitives.c`, `GRHayL/Con2Prim/con2prim_multi_method.c`; fixture: `Unit_Tests/unit_test_con2prim_multi_method_hybrid.c`. |
| Hybrid select and multi-method | Undensitized `cons_undens` | Primitive recovery result | Select dispatches one hybrid key; multi-method owns optional guess and backup loop. Source: `GRHayL/Con2Prim/con2prim_multi_method.c`; tests: `Unit_Tests/unit_test_con2prim_multi_method_hybrid.c`, `Unit_Tests/unit_test_hybrid_failure.c`. |
| Tabulated select and multi-method | Undensitized `cons_undens` | Primitive recovery result | Same selector versus multi-method split as hybrid, but compiled-out HDF5 paths return `ghl_error_used_disabled_hdf5`. Source: `GRHayL/Con2Prim/con2prim_multi_method.c`; tests: `Unit_Tests/unit_test_con2prim_tabulated.c`, `Unit_Tests/unit_test_code_error.c`. |
| `ghl_enforce_primitive_limits_and_compute_u0` | No conservative input | No conservative output | Post-recovery primitive limiter; recomputes EOS fields and `u0`, and reports speed limiting through the caller-provided flag. Source: `GRHayL/Con2Prim/enforce_primitive_limits_and_compute_u0.c`; tests: `Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c`. |
| `ghl_compute_conservs` | Primitives with initialized `eps` and `u0` | Densitized conservative `cons` | Primitive-to-conservative recompute after primitive limits. Source: `GRHayL/Con2Prim/compute_conservs.c`; tests: `Unit_Tests/unit_test_compute_conservs_and_Tmunu.c`. |
| `ghl_compute_conservs_and_Tmunu` | Primitives with initialized `eps` and `u0` | Densitized conservative `cons` and `Tmunu` | Same conservative recompute plus stress-energy tensor. Source: `GRHayL/Con2Prim/compute_conservs_and_Tmunu.c`; tests: `Unit_Tests/unit_test_compute_conservs_and_Tmunu.c`, `Unit_Tests/data_gen/unit_test_data_con2prim_multi_method_hybrid.c`. |

## Test Routes

- Conservative limits: `Unit_Tests/unit_test_apply_conservative_limits.c`.
- Hybrid selected-method recovery and method fixtures:
  `Unit_Tests/unit_test_con2prim_multi_method_hybrid.c`.
- Hybrid multi-method failure and backup behavior:
  `Unit_Tests/unit_test_hybrid_failure.c`.
- Tabulated multi-method recovery, backup accounting, and EOS-table paths:
  `Unit_Tests/unit_test_con2prim_tabulated.c`.
- Primitive post-limits and `u0`:
  `Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c`.
- Conservative recompute and stress-energy:
  `Unit_Tests/unit_test_compute_conservs_and_Tmunu.c`.
