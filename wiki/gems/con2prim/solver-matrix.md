# Con2Prim Solver Matrix

Purpose: source-grounded route for Con2Prim method support, dispatch, build
evidence, public API contracts, and diagnostics. This page covers
`GRHayL/Con2Prim/**`, `GRHayL/include/ghl.h`,
`GRHayL/include/ghl_con2prim.h`, `docs/raw/Con2Prim.dox`, and direct
Con2Prim tests. Repo-local source remains authority.

Support labels below require more than file existence: enum/name mapping,
public declaration, selector dispatch, build-list inclusion, and test evidence
where applicable.

## Method-Support Matrix

| method/key | enum evidence | routine-name mapping | public declaration | selector case | build-list evidence | Doxygen table status | test evidence | status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `None` | `ghl_con2prim_id_None` in `GRHayL/include/ghl.h` | `"None"` in `GRHayL/Con2Prim/get_con2prim_routine_name.c` | No solver declaration; used in `ghl_parameters.backup_routine[3]` | No solver dispatch; multi-method stops backup loop when a backup slot is `None` in `GRHayL/Con2Prim/con2prim_multi_method.c` | No solver build entry | Not listed | Used as no-backup sentinel in `Unit_Tests/unit_test_con2prim_multi_method_hybrid.c`, `Unit_Tests/unit_test_con2prim_tabulated.c`, and `Unit_Tests/unit_test_con2prim_debug.c` | Sentinel only, not solver |
| `Noble2D` | `ghl_con2prim_id_Noble2D` in `GRHayL/include/ghl.h` | `"Noble2D"` in `GRHayL/Con2Prim/get_con2prim_routine_name.c` | `ghl_hybrid_Noble2D` and `ghl_tabulated_Noble2D` in `GRHayL/include/ghl_con2prim.h` | Hybrid and tabulated selector cases in `GRHayL/Con2Prim/con2prim_multi_method.c` | Hybrid: `GRHayL/Con2Prim/Hybrid/Noble/Noble2D/make.code.defn`; tabulated: `GRHayL/Con2Prim/Tabulated/Noble2D/make.code.defn` | Hybrid and tabulated tables list it in `docs/raw/Con2Prim.dox` | Hybrid coverage in `Unit_Tests/unit_test_con2prim_multi_method_hybrid.c`; failure/backup coverage in `Unit_Tests/unit_test_hybrid_failure.c`; tabulated run/data coverage in `Unit_Tests/unit_test_con2prim_tabulated.c` and `.github/run_tests.sh` | Hybrid supported; tabulated supported when HDF5 enabled |
| `Noble1D` | `ghl_con2prim_id_Noble1D` in `GRHayL/include/ghl.h` | `"Noble1D"` in `GRHayL/Con2Prim/get_con2prim_routine_name.c` | `ghl_hybrid_Noble1D` in `GRHayL/include/ghl_con2prim.h` | Hybrid selector case in `GRHayL/Con2Prim/con2prim_multi_method.c`; no tabulated case | `GRHayL/Con2Prim/Hybrid/Noble/Noble1D/make.code.defn` | Hybrid table lists it | `Unit_Tests/unit_test_con2prim_multi_method_hybrid.c` | Hybrid supported |
| `Noble1D_entropy` | `ghl_con2prim_id_Noble1D_entropy` in `GRHayL/include/ghl.h` | `"Noble1D_entropy"` in `GRHayL/Con2Prim/get_con2prim_routine_name.c` | `ghl_hybrid_Noble1D_entropy` in `GRHayL/include/ghl_con2prim.h` | Hybrid selector case in `GRHayL/Con2Prim/con2prim_multi_method.c`; no tabulated case | `GRHayL/Con2Prim/Hybrid/Noble/Noble1D/make.code.defn` | Hybrid table lists it | `Unit_Tests/unit_test_con2prim_multi_method_hybrid.c` | Hybrid supported when entropy path is intended |
| `Noble1D_entropy2` | `ghl_con2prim_id_Noble1D_entropy2` in `GRHayL/include/ghl.h` | `"Noble1D_entropy2"` in `GRHayL/Con2Prim/get_con2prim_routine_name.c` | `ghl_hybrid_Noble1D_entropy2` in `GRHayL/include/ghl_con2prim.h` | No hybrid or tabulated selector case found in `GRHayL/Con2Prim/con2prim_multi_method.c` | Source file `GRHayL/Con2Prim/Hybrid/Noble/Noble1D/hybrid_Noble1D_entropy2.c` exists, but it is absent from `GRHayL/Con2Prim/Hybrid/Noble/Noble1D/make.code.defn` | Not listed | No direct Con2Prim solver test found | Unresolved/source-present: enum, name, declaration, and source exist; selector/build/test evidence missing |
| `Font1D` | `ghl_con2prim_id_Font1D` in `GRHayL/include/ghl.h` | `"Font1D"` in `GRHayL/Con2Prim/get_con2prim_routine_name.c` | `ghl_hybrid_Font1D` in `GRHayL/include/ghl_con2prim.h` | Hybrid selector case in `GRHayL/Con2Prim/con2prim_multi_method.c`; no tabulated case | `GRHayL/Con2Prim/Hybrid/Font1D/make.code.defn` | Hybrid table lists it | `Unit_Tests/unit_test_con2prim_multi_method_hybrid.c`; backup/failure coverage in `Unit_Tests/unit_test_hybrid_failure.c` | Hybrid supported |
| `Palenzuela1D` | `ghl_con2prim_id_Palenzuela1D` in `GRHayL/include/ghl.h` | `"Palenzuela1D"` in `GRHayL/Con2Prim/get_con2prim_routine_name.c` | `ghl_hybrid_Palenzuela1D_energy` and `ghl_tabulated_Palenzuela1D_energy` in `GRHayL/include/ghl_con2prim.h` | Hybrid and tabulated selector cases in `GRHayL/Con2Prim/con2prim_multi_method.c` | Hybrid: `GRHayL/Con2Prim/Hybrid/Palenzuela1D/make.code.defn`; tabulated: `GRHayL/Con2Prim/Tabulated/Palenzuela1D/make.code.defn` | Hybrid and tabulated tables list it | Hybrid coverage in `Unit_Tests/unit_test_con2prim_multi_method_hybrid.c`; tabulated coverage in `Unit_Tests/unit_test_con2prim_tabulated.c` and `.github/run_tests.sh` | Hybrid supported; tabulated supported when HDF5 enabled |
| `Palenzuela1D_entropy` | `ghl_con2prim_id_Palenzuela1D_entropy` in `GRHayL/include/ghl.h` | `"Palenzuela1D_entropy"` in `GRHayL/Con2Prim/get_con2prim_routine_name.c` | `ghl_hybrid_Palenzuela1D_entropy` and `ghl_tabulated_Palenzuela1D_entropy` in `GRHayL/include/ghl_con2prim.h` | Hybrid and tabulated selector cases in `GRHayL/Con2Prim/con2prim_multi_method.c` | Hybrid: `GRHayL/Con2Prim/Hybrid/Palenzuela1D/make.code.defn`; tabulated: `GRHayL/Con2Prim/Tabulated/Palenzuela1D/make.code.defn` | Hybrid and tabulated tables list it | Hybrid coverage in `Unit_Tests/unit_test_con2prim_multi_method_hybrid.c`; tabulated coverage in `Unit_Tests/unit_test_con2prim_tabulated.c` and `.github/run_tests.sh` | Hybrid supported; tabulated supported when HDF5 enabled and entropy path is intended |
| `Newman1D` | `ghl_con2prim_id_Newman1D` in `GRHayL/include/ghl.h` | `"Newman1D"` in `GRHayL/Con2Prim/get_con2prim_routine_name.c` | `ghl_tabulated_Newman1D_energy` in `GRHayL/include/ghl_con2prim.h` | Tabulated selector case in `GRHayL/Con2Prim/con2prim_multi_method.c`; no hybrid case | `GRHayL/Con2Prim/Tabulated/Newman1D/make.code.defn` | Tabulated table lists it | `Unit_Tests/unit_test_con2prim_tabulated.c` uses it with `Palenzuela1D` backup; fixture routes in `.github/run_tests.sh`; debug-only route in `Unit_Tests/unit_test_con2prim_debug.c` | Tabulated supported when HDF5 enabled |
| `Newman1D_entropy` | `ghl_con2prim_id_Newman1D_entropy` in `GRHayL/include/ghl.h` | `"Newman1D_entropy"` in `GRHayL/Con2Prim/get_con2prim_routine_name.c` | `ghl_tabulated_Newman1D_entropy` in `GRHayL/include/ghl_con2prim.h` | Tabulated selector case in `GRHayL/Con2Prim/con2prim_multi_method.c`; no hybrid case | `GRHayL/Con2Prim/Tabulated/Newman1D/make.code.defn` | Tabulated table lists it | `Unit_Tests/unit_test_con2prim_tabulated.c`; fixture routes in `.github/run_tests.sh` | Tabulated supported when HDF5 enabled and entropy path is intended |
| Source-present non-enum Cerda-Duran file | No `ghl_con2prim_id_*` enum found for Cerda-Duran | No `ghl_get_con2prim_routine_name` case found | No public `ghl_*Cerda*` declaration found | No selector dispatch case found | `GRHayL/Con2Prim/Tabulated/con2prim_CerdaDuran3D.cc` exists, but `GRHayL/Con2Prim/Tabulated/make.code.defn` lists only `Newman1D Noble2D Palenzuela1D` subdirs | Not listed | No direct test or CI runner evidence found | Source-present unresolved; not documented as supported |

## `ghl_con2prim_id_t`

`ghl_con2prim_id_t` is declared in `GRHayL/include/ghl.h` and currently
contains `None`, `Noble2D`, `Noble1D`, `Noble1D_entropy`,
`Noble1D_entropy2`, `Font1D`, `Palenzuela1D`, `Palenzuela1D_entropy`,
`Newman1D`, and `Newman1D_entropy`. The enum is a public key space, not proof
that every key is active in every EOS family. Confirm against
`GRHayL/Con2Prim/con2prim_multi_method.c` and the relevant `make.code.defn`
before calling a method supported.

## `ghl_parameters` Con2Prim Fields

Selection fields live in `GRHayL/include/ghl.h` and are initialized by
`GRHayL/GRHayL_Core/initialize_params.c`:

- `main_routine`: primary `ghl_con2prim_id_t` used by the multi-method drivers.
- `backup_routine[3]`: up to three fallback `ghl_con2prim_id_t` values; `None`
  ends backup attempts in `GRHayL/Con2Prim/con2prim_multi_method.c`.
- `evolve_entropy`: caller-side control used by entropy-aware tests and solver
  setup.
- `evolve_temp`: tabulated-temperature evolution flag stored in parameters.
- `calc_prim_guess`: when true, multi-method drivers call
  `ghl_guess_primitives` before selector dispatch.
- `psi6threshold`: affects conservative and primitive limit handling.
- `max_Lorentz_factor` and `inv_sq_max_Lorentz_factor`: speed-limit caps used by
  velocity limiting and related recovery helpers.
- `con2prim_max_iterations` and `con2prim_solver_tolerance`: defaulted in
  `initialize_params.c`; consumed by Noble initialization paths. Current
  Palenzuela and Newman paths use local hard-coded iteration/tolerance values.

Tabulated EOS also carries `root_finding_precision` in `ghl_eos_parameters`
(`GRHayL/include/ghl.h`) and initializes it in `GRHayL/GRHayL_Core/initialize_eos.c`;
tabulated tests override it in `Unit_Tests/unit_test_con2prim_tabulated.c`.

## Routine-Name Mapping

`ghl_get_con2prim_routine_name` is declared in `GRHayL/include/ghl.h` and
implemented in `GRHayL/Con2Prim/get_con2prim_routine_name.c`. It maps enum keys
to display strings and returns `NULL` for unknown keys. This name helper does
not dispatch solvers and does not prove build support.

## Public Select And Multi-Method APIs

Public declarations are in `GRHayL/include/ghl_con2prim.h`:

- `ghl_con2prim_hybrid_select_method`: dispatches one hybrid key. Unknown or
  unsupported keys return `ghl_error_invalid_c2p_key`.
- `ghl_con2prim_tabulated_select_method`: dispatches one tabulated key when HDF5
  is enabled. Under `GHL_DISABLE_HDF5`, it returns `ghl_error_used_disabled_hdf5`.
- `ghl_con2prim_hybrid_multi_method`: optionally calls `ghl_guess_primitives`,
  runs `params->main_routine`, then tries up to three backups.
- `ghl_con2prim_tabulated_multi_method`: same backup structure for tabulated
  methods, with the same HDF5-disabled guard.
- `ghl_con2prim_multi_method`: exported function pointer declaration; confirm
  assignment at integration/configuration sites before relying on it.

Both selector functions receive undensitized conservative quantities by local
parameter naming in `GRHayL/Con2Prim/con2prim_multi_method.c`. Callers that
start from densitized conservatives should route through
`ghl_undensitize_conservatives` first.

## Diagnostics Contracts

`ghl_con2prim_diagnostics` is declared in `GRHayL/include/ghl_con2prim.h` and
initialized by `GRHayL/Con2Prim/initialize_diagnostics.c`.

- `tau_fix`: initialized false; set true by
  `GRHayL/Con2Prim/apply_conservative_limits.c` when conservative energy is
  adjusted.
- `Stilde_fix`: initialized false; set true by
  `GRHayL/Con2Prim/apply_conservative_limits.c` when momentum is rescaled.
- `speed_limited`: initialized false; set by solver finalization paths and by
  `ghl_limit_utilde_and_compute_v` / primitive-limit helpers when speed limiting
  occurs.
- `which_routine`: initialized to `ghl_con2prim_id_None`; set by successful
  solver wrapper paths such as Noble, Font, Palenzuela, and Newman routines.
- `backup[3]`: initialized false; multi-method drivers set slot `n` true before
  trying `params->backup_routine[n]`.
- `n_iter`: declared as number of iterations, but
  `ghl_initialize_diagnostics` does not initialize it. Source shows several
  solvers assign it from local iteration counters (`harm_aux.n_iter`,
  `rparams.n_iters`, or `step`), while `Font1D` does not show an assignment in
  the direct wrapper search. Do not promise stable `n_iter` semantics without
  checking the specific solver source.

## Ground Truth

- `GRHayL/include/ghl.h`
- `GRHayL/include/ghl_con2prim.h`
- `GRHayL/Con2Prim/get_con2prim_routine_name.c`
- `GRHayL/Con2Prim/con2prim_multi_method.c`
- `GRHayL/Con2Prim/**/make.code.defn`
- `GRHayL/Con2Prim/Hybrid/**`
- `GRHayL/Con2Prim/Tabulated/**`
- `docs/raw/Con2Prim.dox`
- `Unit_Tests/unit_test_con2prim_multi_method_hybrid.c`
- `Unit_Tests/unit_test_con2prim_tabulated.c`
- `Unit_Tests/unit_test_hybrid_failure.c`
- `Unit_Tests/unit_test_con2prim_debug.c`
- `.github/run_tests.sh`
