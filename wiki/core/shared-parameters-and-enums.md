# Shared Parameters And Enums

> Shared `ghl.h` parameters, enums, struct routes, `GRHayL_Core/make.code.defn` source mapping, and non-Core `ghl.h` declarations. · Status: confirmed
> Up: [GRHayL Core](index.md)

## Summary

`GRHayL/include/ghl.h` (with `ghl_io.h`, `ghl_debug.h`, and `ghl_metric_helpers.h`) defines GRHayL's shared public enums, parameter structs (`ghl_parameters`, `ghl_eos_parameters`), fluid/metric/stress-energy structs, and normalization macros (`SQR`, `ONE_OVER_SQRT_4PI`, `M_PI`, `MAX_EOS_PARAMS`). Most `ghl.h` declarations are implemented in `GRHayL/GRHayL_Core/`, except the Con2Prim declarations audited below.

## Detail

### `GRHayL_Core` Source Map And Non-Core `ghl.h` Declaration Audit

Every C file listed in `GRHayL/GRHayL_Core/make.code.defn` maps to one defining Core leaf:

| `GRHayL_Core` C file(s) | Defining Core leaf |
| --- | --- |
| `abort_if_error.c`, `info_warn_error.c`, `min_max_clamp.c` | [Errors, IO, Debug, And Utilities](errors-io-debug-utilities.md) (Ev: `GRHayL/GRHayL_Core/make.code.defn::abort_if_error.c`) |
| `compute_ADM_auxiliaries.c`, `enforce_detgij_and_initialize_ADM_metric.c`, `initialize_metric.c` | [Metric And ADM Helpers](metric-adm.md) (Ev: `GRHayL/GRHayL_Core/make.code.defn::initialize_metric.c`) |
| `compute_smallb_and_b2.c`, `compute_TDNmunu.c`, `compute_TUPmunu.c`, `initialize_stress_energy.c`, `return_stress_energy.c` | [Stress-Energy And `smallb`](stress-energy-smallb.md) (Ev: `GRHayL/GRHayL_Core/make.code.defn::compute_TDNmunu.c`) |
| `initialize_conservatives.c`, `initialize_extrinsic_curvature.c`, `initialize_params.c`, `initialize_primitives.c`, `return_conservatives.c`, `return_primitives.c` | [Struct Pack And Unpack](struct-pack-unpack.md) (Ev: `GRHayL/GRHayL_Core/make.code.defn::initialize_params.c`) |
| `initialize_eos.c` | [EOS Dispatch](eos-dispatch.md) (Ev: `GRHayL/GRHayL_Core/make.code.defn::initialize_eos.c`) |
| `limit_v_and_compute_u0.c` | [Velocity And `u0`](velocity-u0.md) (Ev: `GRHayL/GRHayL_Core/make.code.defn::limit_v_and_compute_u0.c`) |

**Non-Core declarations in `ghl.h`:**
- `ghl_get_con2prim_routine_name` is declared in `GRHayL/include/ghl.h` but implemented in `GRHayL/Con2Prim/get_con2prim_routine_name.c`.
- `ghl_compute_SU_Bsq_Ssq_BdotS` is declared in `GRHayL/include/ghl.h` but implemented in `GRHayL/Con2Prim/compute_SU_Bsq_Ssq_BdotS.c` ([Limits And Conversions](../gems/con2prim/limits-and-conversions.md)).
- The opaque `ghl_c2p_nn_model` forward declaration in `ghl.h` belongs to [Neural-Network Primitive Guess](../gems/con2prim/neural-network-primitive-guess.md).

Evidence (descriptive behavior): `ghl.h::ghl_get_con2prim_routine_name`; corroboration: `get_con2prim_routine_name.c` and `compute_SU_Bsq_Ssq_BdotS.c`.

### Public Enums And Shared Struct Routes

| Type | Defining route | Purpose |
| --- | --- | --- |
| `ghl_error_codes_t` / `ghl_error_keys` | [Errors, IO, Debug, And Utilities](errors-io-debug-utilities.md) | Recoverable return codes (`abort_if_error.c` messages) and `ghl_io.h` macro controls. (Ev: `ghl.h::ghl_error_codes_t`) |
| `ghl_con2prim_id_t` | [Solver Matrix](../gems/con2prim/solver-matrix.md) | Solver IDs (`Noble2D`, `Font1D`, `Palenzuela1D`, `Newman1D`, entropy variants, `None`). (Ev: `ghl.h::ghl_con2prim_id_t`) |
| `ghl_eos_t` / `ghl_eos_table_t` | [EOS Initialization And Dispatch](../gems/eos/initialization-and-dispatch.md) | Simple/hybrid/tabulated family enum and table-format enum (`ghl_eos_table_stellarcollapse`). (Ev: `ghl.h::ghl_eos_t`) |
| `ghl_primitive_quantities` / `ghl_conservative_quantities` / `ghl_extrinsic_curvature` | [Struct Pack And Unpack](struct-pack-unpack.md) | Primitive, conservative, and symmetric `K_{ij}` point structs (`ghl.h`). (Ev: `ghl.h::ghl_extrinsic_curvature`) |
| `ghl_metric_quantities` / `ghl_ADM_aux_quantities` | [Metric And ADM Helpers](metric-adm.md) | ADM 3-metric and 4-metric structs (`ghl.h`). (Ev: `ghl.h::ghl_ADM_aux_quantities`) |
| `ghl_stress_energy` | [Stress-Energy And `smallb`](stress-energy-smallb.md) | `T4[4][4]` stress-energy struct (`ghl.h`). (Ev: `ghl.h::ghl_stress_energy`) |

### `ghl_parameters` And `ghl_eos_parameters` Field Groups

- **`ghl_parameters` (`initialize_params.c`):** stores `main_routine`, `backup_routine[3]`, `evolve_entropy`, `evolve_temp`, `calc_prim_guess`, `psi6threshold`, `max_Lorentz_factor` (and derived `inv_sq_max_Lorentz_factor`), `con2prim_max_iterations`, `con2prim_solver_tolerance`, `ppm_*` flattening/shock parameters ([PPM Flow](../gems/reconstruction/ppm-flow.md)), and `Lorenz_damping_factor` ([Gauge RHS](../gems/induction/gauge-rhs.md)). `ghl_initialize_params` returns `void` with no validation, reads all three `backup_routine` entries, and divides `1.0 / max_Lorentz_factor^2` directly; callers must pass finite `max_Lorentz_factor >= 1`.
- **`ghl_eos_parameters` (`initialize_eos.c`):** stores `eos_type`, `table_type`, `clean_sound_speed`, atmosphere/min/max bounds (`rho_*`, `tau_atm`, `press_*`, `Y_e_*`, `T_*`, `eps_*`, `entropy_*`), piecewise-polytrope arrays (`neos`, `rho_ppoly`, `Gamma_ppoly`, `K_ppoly`, `eps_integ_const`, `p_ppoly`, `Gamma_th`), and tabulated table arrays/beta-equilibrium fields/`root_finding_precision`. On success, Core EOS initializers replace the destination with a zero-initialized aggregate populated for the active family (setting unused `Y_e_atm` and `T_atm` to zero for simple/hybrid so constant Atmosphere can copy them).
- **Cross-module edit hazards:** modifying `ghl_con2prim_id_t`, `ghl_error_codes_t`, `ghl_eos_parameters`, `ppm_*`, or `Lorenz_damping_factor` impacts downstream gems, `abort_if_error.c`, and `implementations/GRHayLib/`.

Evidence (descriptive behavior): `initialize_params.c::ghl_initialize_params`; corroboration: `ghl.h::ghl_parameters` and `ghl.h::ghl_eos_parameters`.

## Sources

- [GRHayL/include/ghl.h](../../GRHayL/include/ghl.h) - public enums, structs, macros, and non-Core declarations
- [GRHayL/include/ghl_io.h](../../GRHayL/include/ghl_io.h) - `ghl_error_keys` and IO macros
- [GRHayL/include/ghl_debug.h](../../GRHayL/include/ghl_debug.h) - inline debug print helpers
- [GRHayL/include/ghl_metric_helpers.h](../../GRHayL/include/ghl_metric_helpers.h) - inline 3D/4D metric helpers
- [GRHayL/GRHayL_Core/make.code.defn](../../GRHayL/GRHayL_Core/make.code.defn) - Core C source list
- [GRHayL/GRHayL_Core/initialize_params.c](../../GRHayL/GRHayL_Core/initialize_params.c) - `ghl_initialize_params`
- [GRHayL/GRHayL_Core/initialize_eos.c](../../GRHayL/GRHayL_Core/initialize_eos.c) - Core EOS aggregate initialization
- [GRHayL/Con2Prim/get_con2prim_routine_name.c](../../GRHayL/Con2Prim/get_con2prim_routine_name.c) - Con2Prim-implemented `ghl.h` function
- [GRHayL/Con2Prim/compute_SU_Bsq_Ssq_BdotS.c](../../GRHayL/Con2Prim/compute_SU_Bsq_Ssq_BdotS.c) - Con2Prim-implemented `ghl.h` function

## See Also

- Depends on: [Public API Map](../architecture/public-api-map.md)
- See also: [Struct Pack And Unpack](struct-pack-unpack.md)
- See also: [EOS Dispatch](eos-dispatch.md)
- See also: [GRHayLib Runtime Parameters](../implementations/grhaylib-runtime-parameters.md)
