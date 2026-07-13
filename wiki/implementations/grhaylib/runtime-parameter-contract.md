# GRHayLib Runtime Parameter Contract

Purpose: route GRHayLib runtime lifecycle, global state, EOS setup,
Con2Prim selection, and Cactus parameter plumbing. This page does not replace
GRHayL Core, EOS, Con2Prim, or Reconstruction contracts; it maps how the
Cactus thorn wires those public APIs together.

Primary implementation truth is
[`implementations/GRHayLib/src/initialize_and_shutdown.c`](../../../implementations/GRHayLib/src/initialize_and_shutdown.c),
[`implementations/GRHayLib/param.ccl`](../../../implementations/GRHayLib/param.ccl),
and
[`implementations/GRHayLib/src/GRHayLib.h`](../../../implementations/GRHayLib/src/GRHayLib.h).
Use [Con2Prim solver matrix](../../gems/con2prim/solver-matrix.md) for solver
support truth and
[EOS initialization and dispatch](../../gems/eos/initialization-and-dispatch.md)
for EOS function-pointer and parameter setup truth.

## Runtime Lifecycle

`GRHayLib.h` is the thorn-facing aggregate header. It includes GRHayL public
headers and exports two global pointers:

- `ghl_parameters *ghl_params`
- `ghl_eos_parameters *ghl_eos`

`GRHayLib_paramcheck` validates Cactus parameters before runtime use. It checks
required atmosphere values, min/max sentinels, EOS-specific setup inputs,
Con2Prim/EOS compatibility, and entropy-gated methods.

`GRHayLib_initialize` owns runtime allocation and setup:

1. Allocates `ghl_params` with `malloc(sizeof(ghl_parameters))`.
2. Allocates `ghl_eos` with `malloc(sizeof(ghl_eos_parameters))`.
3. Maps `con2prim_routine` and `con2prim_backup_routines[3]` through
   `parse_C2P_routine_keyword`.
4. Calls `GRHayLib_paramcheck`.
5. Calls `ghl_initialize_params`.
6. Copies PPM Cactus parameters into `ghl_params`.
7. Initializes EOS parameters and EOS-dependent dispatch for `EOS_type`.

Both allocations use unchecked `malloc`, not zero-initializing allocation.
There is no local rollback path if setup later errors. Tabulated initialization
can itself leave partial table state on several failures; see the
[tabulated table contract](../../gems/eos/tabulated-table-contract.md). Treat
successful initialization as a precondition for consumers and termination.

`GRHayLib_terminate` owns shutdown:

1. If `ghl_eos->eos_type == ghl_eos_tabulated`, call
   `ghl_tabulated_free_memory(ghl_eos)` before freeing EOS storage.
2. Free `ghl_eos`.
3. Free `ghl_params`.
4. Set both globals to `NULL`.

No other owner is visible for these two global allocations in GRHayLib source.
The schedule skips `GRHayLib_initialize` when `ID_converter_ILGRMHD` is active,
but still schedules `GRHayLib_terminate`; this repository does not contain that
thorn or prove who initializes these globals on the skipped path. Safe ownership
there needs downstream confirmation before runtime claims.

## `ghl_initialize_params` Map

GRHayLib maps Cactus controls into `ghl_initialize_params` at routing level:

- Main and backup Con2Prim choices: `con2prim_routine` and
  `con2prim_backup_routines[3]`, after `parse_C2P_routine_keyword`.
- Evolution flags: `evolve_entropy` and `evolve_temperature`.
- Primitive guess policy: `calc_primitive_guess`.
- Limit and gauge controls: `Psi6threshold`, `max_Lorentz_factor`, and
  `Lorenz_damping_factor`.

The tabulated NN fallback toggle `enable_backup_nn_primitive_guess` is not
part of `ghl_initialize_params`; GRHayLib passes it to
`ghl_initialize_tabulated_eos` as the `enable_neural_net_c2p` argument.

Core then stores these in `ghl_params` as `main_routine`,
`backup_routine[3]`, `evolve_entropy`, `evolve_temp`, `calc_prim_guess`,
`psi6threshold`, `max_Lorentz_factor`, derived
`inv_sq_max_Lorentz_factor`, and `Lorenz_damping_factor`. Core also sets
default Con2Prim iteration/tolerance and default PPM controls. See
[shared parameters and enums](../../core/shared-parameters-and-enums.md).

## PPM Overrides

After `ghl_initialize_params`, GRHayLib overwrites Core PPM defaults with
Cactus parameter values:

- `ppm_flattening_epsilon`
- `ppm_flattening_omega1`
- `ppm_flattening_omega2`
- `ppm_shock_epsilon`
- `ppm_shock_eta1`
- `ppm_shock_eta2`
- `ppm_shock_k0`

These fields route to [PPM flow](../../gems/reconstruction/ppm-flow.md). Keep
PPM behavior there; this page only records GRHayLib parameter plumbing.

## EOS Branches

`EOS_type = "Simple"`:

- `GRHayLib_paramcheck` requires `Gamma`, `rho_b_atm`, and `P_atm`.
- Optional `rho_b_min`, `rho_b_max`, `P_min`, and `P_max` may use `-1`
  sentinel values.
- `GRHayLib_initialize` sets `ghl_con2prim_multi_method =
  ghl_con2prim_hybrid_multi_method`.
- It calls `ghl_initialize_simple_eos_functions_and_params(rho_b_atm,
  rho_b_min, rho_b_max, P_atm, P_min, P_max, Gamma, ghl_eos)`.
- Core wrapper installs EOS function pointers through
  `ghl_initialize_eos_functions(ghl_eos_simple)`, then calls
  `ghl_initialize_simple_eos`.

`EOS_type = "Hybrid"`:

- `GRHayLib_paramcheck` requires `Gamma_th`, `Gamma_ppoly_in[0]`, and
  `k_ppoly0`; for `neos > 1`, it also checks required `rho_ppoly_in` and
  `Gamma_ppoly_in` entries.
- Hybrid uses `rho_b_atm` and optional `rho_b_min`/`rho_b_max`; pressure,
  energy, entropy, and tau atmosphere/bounds are computed by EOS helpers.
- `GRHayLib_initialize` sets `ghl_con2prim_multi_method =
  ghl_con2prim_hybrid_multi_method`.
- It calls `ghl_initialize_hybrid_eos_functions_and_params(rho_b_atm,
  rho_b_min, rho_b_max, neos, rho_ppoly_in, Gamma_ppoly_in, k_ppoly0,
  Gamma_th, ghl_eos)`.
- Core wrapper installs EOS function pointers through
  `ghl_initialize_eos_functions(ghl_eos_hybrid)`, then calls
  `ghl_initialize_hybrid_eos`.
- Current upstream hybrid initialization can read `rho_ppoly[neos-1]` after
  copying only `neos-1` breakpoints. GRHayLib's CCL bound on `neos` does not
  resolve that upstream source defect.

`EOS_type = "Tabulated"`:

- `GRHayLib_paramcheck` requires `rho_b_atm`, `Y_e_atm`, and `T_atm`, and
  accepts optional `rho_b_min`/`rho_b_max`, `Y_e_min`/`Y_e_max`, and
  `T_min`/`T_max`.
- `GRHayLib_initialize` requires nonempty `EOS_tablepath`.
- It sets `ghl_con2prim_multi_method =
  ghl_con2prim_tabulated_multi_method`.
- It calls `ghl_initialize_eos_functions(ghl_eos_tabulated)` directly.
- It then calls `ghl_initialize_tabulated_eos(EOS_tablepath,
  parse_eos_table_type_keyword(eos_table_type),
  eos_table_clean_sound_speed, enable_backup_nn_primitive_guess, rho_b_atm,
  rho_b_min, rho_b_max, Y_e_atm, Y_e_min, Y_e_max, T_atm, T_min, T_max,
  ghl_eos)`.
- `parse_eos_table_type_keyword` maps `eos_table_type = "stellarcollapse"` to
  `ghl_eos_table_stellarcollapse`; unknown values return
  `ghl_eos_table_unknown`.
- `eos_table_clean_sound_speed` passes the clean-sound-speed flag into
  `ghl_initialize_tabulated_eos`.
- `enable_backup_nn_primitive_guess` passes the downstream neural-network
  primitive-guess flag into `ghl_initialize_tabulated_eos` as
  `enable_neural_net_c2p`.
- Tabulated min/max defaults are table-bound aware: requested lower bounds
  below table limits are replaced by table lower bounds, and missing or
  oversized upper bounds are replaced by table upper bounds.

`EOS_type = ""` errors as unset. Other strings error as unsupported by
`param.ccl`.

## Required And Optional Limits

Required atmosphere values:

- All EOS families require `rho_b_atm`.
- Simple additionally requires `P_atm`.
- Tabulated additionally requires `Y_e_atm` and `T_atm`.

Optional limits:

- Simple: `rho_b_min`, `rho_b_max`, `P_min`, and `P_max` may stay `-1`; Core
  defaults density and pressure floors to `0.0` and ceilings to `1e300`.
- Hybrid: `rho_b_min` and `rho_b_max` may stay `-1`; Core defaults density
  floor to `0.0` and ceiling to `1e300`. Pressure, epsilon, entropy, and tau
  bounds are computed from hybrid EOS data.
- Tabulated: `rho_b_min`, `rho_b_max`, `Y_e_min`, `Y_e_max`, `T_min`, and
  `T_max` are clamped or defaulted against table bounds after table read.
  Pressure, epsilon, and entropy bounds come from table metadata.

## Con2Prim Routing

`parse_C2P_routine_keyword` maps parser inputs to `ghl_con2prim_id_t` values
for `con2prim_routine` and `con2prim_backup_routines[3]`; `param.ccl` exposes
those inputs differently for main and backup parameters:

- `None` is exposed only for `con2prim_backup_routines[3]`; the parser still
  maps it to `ghl_con2prim_id_None`.
- `Noble1D_entropy2` is a parser case, but it is commented out in both local
  `param.ccl` keyword lists.
- Exposed method strings include `Noble2D`, `Noble1D`, `Noble1D_entropy`,
  `Font1D`, `Palenzuela1D`, `Palenzuela1D_entropy`, `Newman1D`, and
  `Newman1D_entropy`.
- Unknown strings return `-100`.

`GRHayLib_paramcheck` error formatting has a separate backup-index defect: each
loop identifies `con2prim_backup_routines[i]` in the label but prints
`con2prim_backup_routines[0]` as the value. For failures at index 1 or 2, the
message can therefore report the wrong routine. Use index and parser evidence,
not the printed value alone, when diagnosing a rejected backup.

Entropy methods require `evolve_entropy`. `GRHayLib_paramcheck` rejects
`Noble1D_entropy`, `Palenzuela1D_entropy`, and `Newman1D_entropy` as main or
backup choices when `evolve_entropy` is false.

Parser/parameter keyword existence does not imply supported solver. Compare
any GRHayLib keyword against [Con2Prim solver matrix](../../gems/con2prim/solver-matrix.md),
selector dispatch, build-list evidence, and tests before claiming support.
Recovery order and backup behavior route to
[Con2Prim recovery flow](../../gems/con2prim/recovery-flow.md).

### Tabulated Neural-Network Primitive Guess

GRHayLib exposes `enable_backup_nn_primitive_guess` as a Cactus parameter. In
the tabulated EOS branch, `GRHayLib_initialize` passes it into
`ghl_initialize_tabulated_eos` as the core `enable_neural_net_c2p` field on
`ghl_eos_parameters`. Core tabulated EOS initialization stores the flag and, if
enabled, loads the embedded NN model from the EOS HDF5 table.

This route applies only to tabulated EOS initialization and the tabulated
Con2Prim fallback-guess path. It is not a separate `ghl_con2prim_id_t` solver,
does not alter simple or hybrid EOS setup, and does not imply simple/hybrid NN
support. The tabulated multi-method still runs the selected main and backup
Con2Prim routines; after failure, `enable_neural_net_c2p` allows retry with
`ghl_c2p_nn_guess_primitives`.

The focused Con2Prim route is
[neural-network primitive guess](../../gems/con2prim/neural-network-primitive-guess.md).
GRHayLib direct-compile build routing also includes
`Con2Prim/Tabulated/neural_network_guess` in
[`implementations/GRHayLib/src/make.code.defn`](../../../implementations/GRHayLib/src/make.code.defn).

## Drift Evidence

Record these source facts as drift evidence only; this page does not propose
source changes:

- `param.ccl` comments out `Noble1D_entropy2` in both `con2prim_routine` and
  `con2prim_backup_routines[3]`.
- `parse_C2P_routine_keyword` still has a parser case for
  `Noble1D_entropy2`.
- `GRHayLib_paramcheck` contains compatibility checks for
  `Newman1D_energy`, while exposed Cactus keywords include `Newman1D` and
  `Newman1D_entropy`. Consequence: `EOS_type = "Simple"` or `"Hybrid"` with
  `con2prim_routine = "Newman1D"` passes paramcheck, then fails at runtime
  on every point through the hybrid selector's
  `ghl_error_invalid_c2p_key` path (no hybrid `Newman1D` case in
  `GRHayL/Con2Prim/con2prim_multi_method.c`). Needs a maintainer decision:
  fix the paramcheck strings or restrict the `param.ccl` keyword lists.
- `schedule.ccl` conditionally skips initialization for
  `ID_converter_ILGRMHD` but always schedules termination. Local files do not
  establish alternate allocation ownership or a safe terminate precondition.
- Simple/hybrid upstream EOS initialization does not set `Y_e_atm` or `T_atm`,
  while the built constant-atmosphere routine copies both fields. GRHayLib does
  not seed the allocations before EOS setup, so those two fields remain
  indeterminate on its simple/hybrid path before atmosphere reset.
- `GRHayLib_terminate` calls tabulated cleanup only after successful state is
  assumed; neither global-pointer null checks nor partial-initialization guards
  are present.

## Local Ground Truth

- [`wiki/public-api-map.md`](../../public-api-map.md)
- [`wiki/core/shared-parameters-and-enums.md`](../../core/shared-parameters-and-enums.md)
- [`wiki/core/eos-dispatch-contract.md`](../../core/eos-dispatch-contract.md)
- [`wiki/gems/con2prim/solver-matrix.md`](../../gems/con2prim/solver-matrix.md)
- [`wiki/gems/con2prim/recovery-flow.md`](../../gems/con2prim/recovery-flow.md)
- [`wiki/gems/eos/initialization-and-dispatch.md`](../../gems/eos/initialization-and-dispatch.md)
- [`wiki/gems/reconstruction/ppm-flow.md`](../../gems/reconstruction/ppm-flow.md)
- [`implementations/GRHayLib/param.ccl`](../../../implementations/GRHayLib/param.ccl)
- [`implementations/GRHayLib/src/GRHayLib.h`](../../../implementations/GRHayLib/src/GRHayLib.h)
- [`implementations/GRHayLib/src/initialize_and_shutdown.c`](../../../implementations/GRHayLib/src/initialize_and_shutdown.c)
- [`implementations/GRHayLib/src/make.code.defn`](../../../implementations/GRHayLib/src/make.code.defn)
- [`implementations/GRHayLib/doc/documentation.tex`](../../../implementations/GRHayLib/doc/documentation.tex)
- [`GRHayL/include/ghl.h`](../../../GRHayL/include/ghl.h)
- [`GRHayL/include/ghl_con2prim.h`](../../../GRHayL/include/ghl_con2prim.h)
- [`GRHayL/include/ghl_eos_functions.h`](../../../GRHayL/include/ghl_eos_functions.h)
