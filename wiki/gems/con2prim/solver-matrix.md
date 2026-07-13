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

| method/key | enum/name/declaration/definition | selector/build | Doxygen/test evidence | GRHayLib seam | status |
| --- | --- | --- | --- | --- | --- |
| `None` | Enum and `"None"` name exist; no solver declaration or definition. | No dispatch/build entry; first `None` backup stops retries. | Not in Doxygen solver tables; tests use it as no-backup sentinel. | Allowed only for backup keyword; parser maps it. | Sentinel only, not solver. |
| `Noble2D` | Enum/name plus hybrid and tabulated declarations and definitions exist. | Both selector cases and both manifest entries exist. | Hybrid Doxygen table lists it. Tabulated Doxygen table omits it despite direct tabulated test coverage. Hybrid selected-method, hybrid failure/backup-attempt, and tabulated tests exercise it. | Main/backup keywords and parser accept it; family checks allow it for simple, hybrid, and tabulated EOS. | Hybrid/simple supported; tabulated supported when HDF5 enabled. Doxygen tabulated omission is documentation drift. |
| `Noble1D` | Enum/name plus hybrid declaration and definition exist. | Hybrid selector/build only. | Hybrid Doxygen table and hybrid selected-method test list it. | Main/backup keywords and parser accept it; tabulated parameter check rejects it. | Hybrid/simple supported. |
| `Noble1D_entropy` | Enum/name plus hybrid declaration and definition exist. | Hybrid selector/build only. | Hybrid Doxygen table and hybrid selected-method test list it. | Main/backup keywords and parser accept it; tabulated check rejects it and GRHayLib requires `evolve_entropy`. | Hybrid/simple supported with valid entropy conservative input; Core dispatch does not gate on `params->evolve_entropy`. |
| `Noble1D_entropy2` | Enum/name/declaration and source definition exist. | No selector case; definition source and sibling `func_rho2.c` are absent from manifest. | Absent from Doxygen and direct tests. | Keyword entries are commented out, but parser still maps the string; family/entropy checks do not classify it. | Source-present unresolved; not supported by current build/dispatch. |
| `Font1D` | Enum/name plus hybrid declaration and definition exist. | Hybrid selector/build only. | Hybrid Doxygen table and hybrid selected-method test list it. `unit_test_hybrid_failure.c` does **not** use Font; it retries Noble2D. | Main/backup keywords and parser accept it; simple and tabulated checks reject it. | Hybrid supported; raw Doxygen says simple excludes it. |
| `Palenzuela1D` | Enum/name plus hybrid/tabulated energy declarations and definitions exist. | Both selector/build paths exist. | Both Doxygen tables and hybrid/tabulated tests list it. | Main/backup keywords and parser accept it for all EOS families. | Hybrid/simple supported; tabulated supported when HDF5 enabled. |
| `Palenzuela1D_entropy` | Enum/name plus hybrid/tabulated entropy declarations and definitions exist. | Both selector/build paths exist. | Both Doxygen tables and hybrid/tabulated tests list it. | Main/backup keywords and parser accept it; GRHayLib requires `evolve_entropy`. | Hybrid/simple supported; tabulated supported when HDF5 is enabled. Both require valid entropy conservative input; Core dispatch does not gate on `params->evolve_entropy`. |
| `Newman1D` | Enum/name plus tabulated energy declaration and definition exist. | Tabulated selector/build only. | Tabulated Doxygen table and tabulated tests list it; debug binary uses it. | Main/backup keywords and parser accept `Newman1D`, but simple/hybrid checks test nonexistent keyword `Newman1D_energy`; incompatible selection can pass parameter validation and later fail hybrid selector dispatch. | Tabulated supported when HDF5 enabled; downstream validation drift. |
| `Newman1D_entropy` | Enum/name plus tabulated entropy declaration and definition exist. | Tabulated selector/build only. | Tabulated Doxygen table and tabulated tests list it. | Main/backup keywords and parser accept it; simple/hybrid checks reject it and GRHayLib requires `evolve_entropy`. | Tabulated supported when HDF5 is enabled with valid entropy conservative input; Core dispatch does not gate on `params->evolve_entropy`. |
| Source-present Cerda-Duran path | No enum/name/public declaration; source contains old-style definitions using identifiers that do not match current public types. | No selector; `con2prim_CerdaDuran3D.cc` is absent from manifest. | Absent from Doxygen/tests. | Cerda-Duran keyword entries are commented out and parser has no case. | Source-present unresolved; no configured compile/support evidence. |

All matrix seams come from `GRHayL/include/ghl.h`,
`GRHayL/include/ghl_con2prim.h`,
`GRHayL/Con2Prim/get_con2prim_routine_name.c`,
`GRHayL/Con2Prim/con2prim_multi_method.c`, recursive Con2Prim
`make.code.defn` files, `docs/raw/Con2Prim.dox`, direct tests,
`implementations/GRHayLib/param.ccl`, and
`implementations/GRHayLib/src/initialize_and_shutdown.c`. GRHayLib mismatches
describe downstream evidence; this page does not assign maintainer intent.

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
  assignment before relying on it. `ghl_initialize_eos_functions` assigns the
  hybrid driver for simple/hybrid EOS and the tabulated driver for tabulated
  EOS in HDF5-enabled builds. GRHayLib also assigns it explicitly. This global
  pointer is process-wide state, not a member of `ghl_eos_parameters`.

Both selector functions receive undensitized conservative quantities by local
parameter naming in `GRHayL/Con2Prim/con2prim_multi_method.c`. Callers that
start from densitized conservatives should route through
`ghl_undensitize_conservatives` first.

## Diagnostics Contracts

`ghl_con2prim_diagnostics` is declared in `GRHayL/include/ghl_con2prim.h` and
initialized by `GRHayL/Con2Prim/initialize_diagnostics.c`.

- `tau_fix`: initialized false; set true only for the later magnetic-energy
  or high-`psi6` corrections in
  `GRHayL/Con2Prim/apply_conservative_limits.c`. The initial unconditional
  `fmax` floor can raise `cons->tau` without setting this flag.
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
