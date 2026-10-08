# Con2Prim Recovery Flow

This page summarizes runtime recovery order and conservative-variable data
boundaries for `GRHayL/Con2Prim/`. It is a route into source, Doxygen, and
tests; repo-local source remains authority.

## Ordered Flow

1. Initialize `ghl_con2prim_diagnostics` with
   `ghl_initialize_diagnostics`. It clears `tau_fix`, `Stilde_fix`,
   `speed_limited`, all three `backup` slots, and sets `which_routine` to
   `ghl_con2prim_id_None`. It also clears `nn_guess_used` and initializes
   `n_iter` to zero. The diagnostics struct and public helper list are
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
   successful solver updates `diagnostics->which_routine` and its applicable
   `diagnostics->n_iter`. Attempt flags are sticky over the logical recovery:
   `speed_limited`, `backup[]`, and `nn_guess_used` record whether the event
   occurred in any attempt, including failed attempts.
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
guess, sets `diagnostics->nn_guess_used = true`, and retries the main selector.
If later tabulated backup attempts fail, the driver reuses the stored
`prims_guess_nn` and calls the same backup selector again. `nn_guess_used`
records that the main NN retry was attempted; it is not set again per backup
and does not prove the NN retry succeeded.

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
- Diagnostics lifecycle is split. Callers initialize diagnostics once per
  logical recovery; conservative limits accumulate `tau_fix` and `Stilde_fix`;
  multi-method records attempted backup slots and NN retry and accumulates
  solver `speed_limited` results; post-limit speed limiting also accumulates it.
  A successful direct solver alone reports only its current call. The successful solver
  owns `which_routine` and its applicable `n_iter`; when no solver succeeds,
  `which_routine` remains `None` and `n_iter` is unspecified.
  Source: `GRHayL/include/ghl_con2prim.h`,
  `GRHayL/Con2Prim/initialize_diagnostics.c`,
  `GRHayL/Con2Prim/apply_conservative_limits.c`,
  `GRHayL/Con2Prim/enforce_primitive_limits_and_compute_u0.c`,
  solver files under `GRHayL/Con2Prim/Hybrid/` and
  `GRHayL/Con2Prim/Tabulated/`.

## Candidate Assessment

This is the test step of first-order flux correction (FOFC): estimate the next
step, flag a cell whose new state requires a floor or fails the
conserved-to-primitive inversion, and recompute the fluxes of flagged cells with
first-order reconstruction. The method is that of Lemaster & Stone (2009) as used
in AthenaK by Fields et al. (2025, Sec. 3.3); see
[Ground Truth References](#ground-truth-references). The relaxed discrete maximum
principle of that paper (its Eq. 12) is not applied here.

`ghl_assess_candidate_state` (`GRHayL/Con2Prim/assess_candidate_state.c`) runs the
flow above on copies of a densitized candidate divided by the cell's
coordinate volume, with the driver chosen by `eos->eos_type`, and sets a flag when
the candidate would need a repair: its `D`, `tau`, or `S_i` is not finite (or its
entropy, when `evolve_entropy` is set, which an energy-based recovery would
otherwise discard) or too large to form the closure bound, recovery fails
numerically, `tau_fix`, `Stilde_fix`, or `speed_limited` is set, Font1D succeeds,
or closure fails. Closure rebuilds the conservatives from the limited primitives
and compares them with the original candidate, which is how primitive limits and
solver clamps larger than the closure bound are seen; a smaller one is flagged
only if a named diagnostic also records it, and a primitive floor with no record
goes unflagged. The evolved entropy is not compared with the energy; with an
entropy-based solver it drives the recovery, so a disagreement appears as a
closure failure. `ghl_apply_conservative_limits` is skipped for tabulated EOS,
because a negative `tau_atm` makes its momentum rescaling take the square root
of a negative number. Configuration errors that are the final result of the
recovery are returned and leave the flag unchanged; an error that a later
backup recovers from is not reported, and a solver the sequence never calls is
not checked, so the caller validates the whole configured solver list at
initialization.

Sources: `GRHayL/Con2Prim/assess_candidate_state.c`,
`GRHayL/Con2Prim/apply_conservative_limits.c`,
`GRHayL/Con2Prim/con2prim_multi_method.c`.
Test: `Unit_Tests/unit_test_flux_correction.c`.

## HDF5-Disabled Tabulated Behavior

When built with `GHL_DISABLE_HDF5`, tabulated select and tabulated
multi-method return `ghl_error_used_disabled_hdf5` before solver dispatch or
backup behavior. The public direct tabulated solver symbols remain
link-visible as non-mutating stubs returning the same error. The tabulated branch in
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
| `ghl_guess_primitives` | Undensitized `cons_undens` by source contract | Primitive guess only | Multi-method drivers and the hybrid selected-method test pass undensitized data. Source: `GRHayL/Con2Prim/guess_primitives.c`, `GRHayL/Con2Prim/con2prim_multi_method.c`; fixture: `Unit_Tests/unit_test_con2prim_multi_method_hybrid.c`. |
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

## Ground Truth References

- J. Fields, H. Zhu, D. Radice, J. M. Stone, W. Cook, S. Bernuzzi, and B. Daszuta,
  "Performance-Portable Binary Neutron Star Mergers with AthenaK", ApJS 276, 35
  (2025): [arXiv:2409.10384](https://arxiv.org/abs/2409.10384). Sec. 3.3 gives
  the FOFC procedure and the relaxed discrete maximum principle (Eqs. 10 to 12).
- M. N. Lemaster and J. M. Stone, "Dissipation and Heating in Supersonic
  Hydrodynamic and MHD Turbulence", ApJ 691, 1092 (2009):
  [doi:10.1088/0004-637X/691/2/1092](https://doi.org/10.1088/0004-637X/691/2/1092),
  the FOFC reference given by Fields et al.
- J. M. Stone et al., "AthenaK: A Performance-Portable Version of the Athena++ AMR
  Framework": [arXiv:2409.16053](https://arxiv.org/abs/2409.16053). The AthenaK
  FOFC implementation is `MHD::FOFC` in `src/mhd/mhd_fofc.cpp` of the
  [AthenaK repository](https://github.com/IAS-Astrophysics/athenak). No AthenaK
  source code is included in GRHayL.
