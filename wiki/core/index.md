# GRHayL Core

> Router for `GRHayL_Core` shared structs, packing helpers, ADM metrics, velocity limiting, stress-energy, EOS dispatch, errors/IO, and Core tests.

| Page | Go here when... |
| --- | --- |
| [Shared Parameters And Enums](shared-parameters-and-enums.md) | You need `ghl_parameters`, `ghl_eos_parameters`, public enums, `GRHayL_Core/make.code.defn` source routing, or the non-Core `ghl.h` declaration audit. |
| [Struct Pack And Unpack](struct-pack-unpack.md) | You need `ghl_initialize_*` / `ghl_return_*` copy routines, `ghl_initialize_params` Con2Prim/PPM defaults, or representation-neutral conservative storage. |
| [Metric And ADM Helpers](metric-adm.md) | You need `ghl_initialize_metric`, `ghl_compute_ADM_auxiliaries`, `ghl_enforce_detgtij_and_initialize_ADM_metric`, or inline `ghl_metric_helpers.h` routines. |
| [Velocity And `u0`](velocity-u0.md) | You need `ghl_limit_v_and_compute_u0`, Lorentz-factor velocity rescaling, 64-retry rounding behavior, or `ghl_error_u0_singular` bounds. |
| [Stress-Energy And `smallb`](stress-energy-smallb.md) | You need `ghl_compute_smallb_and_b2`, `ghl_compute_TDNmunu`, `ghl_compute_TUPmunu`, upper-triangle write semantics, or stress-energy packing helpers. |
| [EOS Dispatch](eos-dispatch.md) | You need `initialize_eos.c` wrappers, process-wide EOS function-pointer assignment, initialization ordering, or `GHL_DISABLE_HDF5` dispatch behavior. |
| [Errors, IO, Debug, And Utilities](errors-io-debug-utilities.md) | You need `ghl_abort_if_error`, `ghl_Warn_Error`, `ghl_debug.h` print helpers, or `ghl_imin`/`ghl_imax`/`ghl_iclamp`/`ghl_clamp`. |
| [Core Tests And Reference Data](tests-and-reference-data.md) | You need `unit_test_grhayl_core_test_suite.c`, Core reference-data generation, perturbation/randomizer helpers, or Core coverage gaps. |

Back to [AGENTS.md](../../AGENTS.md).
