# Con2Prim Solver Matrix

> Method-by-method support across `ghl_con2prim_id_t`, selector dispatch, manifests, tests, GRHayLib, and `ghl_con2prim_diagnostics`. · Status: confirmed
> Up: [Con2Prim](index.md)

## Summary

A Con2Prim method is supported only when enum/name mapping, public declaration, selector dispatch, `make.code.defn` inclusion, and tests agree. `None` is a backup-termination sentinel, `Noble1D_entropy` and `Noble1D_entropy2` require one-piece constant-Gamma EOS (`Gamma_th == Gamma_ppoly[0]`), tabulated methods require HDF5, and `con2prim_CerdaDuran3D.cc` is an unbuilt archival reference.

## Detail

### Method Support Across EOS Families And GRHayLib

| Method / key | Hybrid / Simple | Tabulated (HDF5) | Dispatch, tests & GRHayLib notes |
| --- | --- | --- | --- |
| `None` | Sentinel only | Sentinel only | First `None` in `backup_routine[3]` stops retries; GRHayLib permits `None` for backup keywords only. (Ev: `GRHayL/Con2Prim/con2prim_multi_method.c::ghl_con2prim_id_None`) |
| `Noble2D` | Supported | Supported | Both selectors, manifests, Doxygen, hybrid/tabulated tests, and GRHayLib accept it. (Ev: `GRHayL/include/ghl_con2prim.h::ghl_hybrid_Noble2D`) |
| `Noble1D` | Supported | Unsupported | Hybrid selector/manifest/tests only; GRHayLib rejects it for tabulated EOS. (Ev: `GRHayL/include/ghl_con2prim.h::ghl_hybrid_Noble1D`) |
| `Noble1D_entropy` | One-piece constant-`Gamma` only | Unsupported | Requires `neos == 1` and `Gamma_th == Gamma_ppoly[0]` (else `ghl_error_invalid_eos_type`) and valid conserved entropy (unchecked; bad entropy fails in recovery, e.g. `ghl_error_neg_pressure`); Core dispatch does not check `params->evolve_entropy` (GRHayLib does). (Ev: `GRHayL/Con2Prim/Hybrid/Noble/Noble1D/hybrid_Noble1D_entropy.c`) |
| `Noble1D_entropy2` | One-piece constant-`Gamma` only | Unsupported | Direct density residual (`func_rho2.c`) returning positive unclamped root (post-limiter enforces density bounds); covered in `unit_test_hybrid_failure.c`; same GRHayLib restrictions. (Ev: `GRHayL/Con2Prim/Hybrid/Noble/Noble1D/hybrid_Noble1D_entropy2.c`) |
| `Font1D` | Hybrid supported | Unsupported | Cold `P`/`eps` use EOS-clamped density while `P/rho` enthalpy closure uses recovered conservation-driven density; `Con2Prim.dox` and GRHayLib exclude simple/tabulated EOS. (Ev: `GRHayL/Con2Prim/Hybrid/Font1D/hybrid_Font1D.c`) |
| `Palenzuela1D` | Supported | Supported | Both selectors, manifests, tests, and GRHayLib accept it for every EOS family. (Ev: `GRHayL/include/ghl_con2prim.h::ghl_tabulated_Palenzuela1D_energy`) |
| `Palenzuela1D_entropy` | Supported | Supported | Hybrid requires conserved entropy; tabulated requires it only when `params->evolve_temp` is true; Core does not check `evolve_entropy` (GRHayLib does). (Ev: `GRHayL/Con2Prim/Tabulated/Palenzuela1D/tabulated_Palenzuela1D_entropy.c`) |
| `Newman1D` / `Newman1D_entropy` | Unsupported | Supported | Tabulated selector/manifest/tests only (`unit_test_con2prim_debug.c` uses `Newman1D`); GRHayLib rejects simple/hybrid and requires `evolve_entropy` for entropy variant. (Ev: `GRHayL/include/ghl_con2prim.h::ghl_tabulated_Newman1D_energy`) |
| Archival `con2prim_CerdaDuran3D.cc` | Unsupported | Unsupported | Unbuilt archival C++ file excluded by `Doxyfile`, manifests, selector, tests, and commented GRHayLib keywords. (Ev: `GRHayL/Con2Prim/Tabulated/con2prim_CerdaDuran3D.cc`) |

### Parameters, Selectors, And Diagnostics

- **`ghl_parameters` controls:** `main_routine`, `backup_routine[3]`, `evolve_entropy`, `evolve_temp` (governs tabulated Palenzuela only, not Newman/Noble), `calc_prim_guess`, `psi6threshold`, `max_Lorentz_factor`/`inv_sq_max_Lorentz_factor`, and `con2prim_max_iterations`/`con2prim_solver_tolerance` (read by Noble; Font, Palenzuela, and Newman use local hard-coded limits). Tabulated EOS separately stores `root_finding_precision` in `ghl_eos_parameters` (`initialize_eos.c`, overridden in `unit_test_con2prim_tabulated.c`).
- **Routine names & dispatch:** `ghl_get_con2prim_routine_name` (`get_con2prim_routine_name.c`) maps enum values to strings (`NULL` on unknown) without dispatching. Single-shot selectors (`ghl_con2prim_{hybrid,tabulated}_select_method`) take **undensitized** `cons_undens` and return `ghl_error_invalid_c2p_key` (or `ghl_error_used_disabled_hdf5` under `GHL_DISABLE_HDF5`). Multi-method drivers (`ghl_con2prim_{hybrid,tabulated}_multi_method`) wrap optional guessing and backup retries; process-wide pointer `ghl_con2prim_multi_method` is assigned by `ghl_initialize_eos_functions` or GRHayLib.
- **`ghl_con2prim_diagnostics` (`initialize_diagnostics.c`):** `tau_fix` and `Stilde_fix` start `false` and accumulate `true` in `apply_conservative_limits.c`; `speed_limited` is written per direct call and OR-accumulated across multi-method attempts; `backup[3]` marks attempted backup slots; `which_routine` starts `ghl_con2prim_id_None` and is written only on solver success; `n_iter` starts `0` with solver-specific semantics ([Internal Numerics And Root-Finding](internal-numerics-and-root-finding.md)).

Evidence (descriptive behavior): `GRHayL/Con2Prim/con2prim_multi_method.c::ghl_con2prim_hybrid_select_method`; corroboration: `GRHayL/Con2Prim/get_con2prim_routine_name.c::ghl_get_con2prim_routine_name` and `GRHayL/Con2Prim/initialize_diagnostics.c::ghl_initialize_diagnostics`.

## Sources

- [GRHayL/include/ghl.h](../../../GRHayL/include/ghl.h) - `ghl_con2prim_id_t`, `ghl_parameters`, `ghl_eos_parameters.root_finding_precision`
- [GRHayL/include/ghl_con2prim.h](../../../GRHayL/include/ghl_con2prim.h) - solver, selector, multi-method, and diagnostics declarations
- [GRHayL/Con2Prim/con2prim_multi_method.c](../../../GRHayL/Con2Prim/con2prim_multi_method.c) - selectors and multi-method dispatch
- [GRHayL/Con2Prim/get_con2prim_routine_name.c](../../../GRHayL/Con2Prim/get_con2prim_routine_name.c) - enum-to-string helper
- [GRHayL/Con2Prim/initialize_diagnostics.c](../../../GRHayL/Con2Prim/initialize_diagnostics.c) - diagnostics initialization
- [GRHayL/Con2Prim/Hybrid/Noble/Noble1D/hybrid_Noble1D_entropy.c](../../../GRHayL/Con2Prim/Hybrid/Noble/Noble1D/hybrid_Noble1D_entropy.c) - constant-Gamma restriction
- [GRHayL/Con2Prim/Hybrid/Noble/Noble1D/hybrid_Noble1D_entropy2.c](../../../GRHayL/Con2Prim/Hybrid/Noble/Noble1D/hybrid_Noble1D_entropy2.c) - density-residual entropy solver
- [GRHayL/Con2Prim/Hybrid/Font1D/hybrid_Font1D.c](../../../GRHayL/Con2Prim/Hybrid/Font1D/hybrid_Font1D.c) - hybrid Font1D solver
- [GRHayL/Con2Prim/Tabulated/Palenzuela1D/tabulated_Palenzuela1D_entropy.c](../../../GRHayL/Con2Prim/Tabulated/Palenzuela1D/tabulated_Palenzuela1D_entropy.c) - `evolve_temp` callback behavior
- [GRHayL/Con2Prim/Tabulated/con2prim_CerdaDuran3D.cc](../../../GRHayL/Con2Prim/Tabulated/con2prim_CerdaDuran3D.cc) - archival unbuilt reference
- [docs/raw/Con2Prim.dox](../../../docs/raw/Con2Prim.dox) - solver tables

## See Also

- Depends on: [Recovery Flow](recovery-flow.md)
- Depends on: [Internal Numerics And Root-Finding](internal-numerics-and-root-finding.md)
- See also: [GRHayLib Runtime Parameters](../../implementations/grhaylib-runtime-parameters.md)
- Validated by: [Con2Prim Tests And Reference Data](tests-and-reference-data.md)
