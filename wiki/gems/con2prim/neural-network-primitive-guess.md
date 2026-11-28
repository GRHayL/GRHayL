# Neural-Network Primitive Guess

> Tabulated Con2Prim neural-network retry API, `{q, r, s, t}` inference and overflow fallbacks, HDF5 model layout (`grhayl_nn_c2p`), memory cleanup, and GRHayLib parameter routing. · Status: confirmed
> Up: [Con2Prim](index.md)

## Summary

GRHayL's neural-network Con2Prim feature predicts a scalar `x` to build an alternate tabulated primitive guess when normal tabulated recovery fails. It is not a standalone `ghl_con2prim_id_t` solver and applies only to tabulated EOS (`eos->enable_neural_net_c2p`), with `GHL_NN_C2P_API_VERSION = 4u` fixing the `{q, r, s, t}` input layout (`q_idx = 0`, `s_idx = 2`). No training or model-accuracy source exists in this repository.

## Detail

### Public API, Inference Arithmetic, And Fallbacks

`GRHayL/include/ghl_con2prim.h` defines `GHL_NN_C2P_API_VERSION` (`4u`, informational—loaders validate fields without negotiating an on-disk version dataset), `ghl_nn_c2p_input_t` (`float q, r, s, t`), `ghl_nn_c2p_guess_t` (`float x`), `ghl_c2p_nn_model`, and `ghl_c2p_nn_{guess,guess_primitives,validate_model,load_hdf5,load_from_eos_hdf5,free}`. Shared state (`enable_neural_net_c2p`, `c2p_nn`) and `ghl_error_nn_c2p_*` codes live in `ghl.h`.

- **Auxiliaries (`tabulated_primitive_guess_helpers.c`):** forms `q = tau / D`, `r = S_squared / D^2`, `s = B_squared / D`, `t = BdotS / D^(3/2)` (`D = cons_undens->rho`).
- **Scalar prediction (`c2p_nn_guess_x.c`, `ghl_c2p_nn.h`):** expects an already-validated model. If `model == NULL`, `q` or `s` is nonfinite, or `x_lo = 1 + q - s` or `width = 1 + q` is nonfinite, it returns zero-initialized `{0}`. Otherwise it floors `width` by `model->dx_eps`, computes fallback `x = x_lo + 0.5 * width` (returning `{0}` if midpoint overflows), applies input transforms (kind `0` identity, kind `1` `log10` after `x_eps` floor), clips to scaling intervals, evaluates hard-tanh hidden layers and sigmoid output clamped to `[y_eps, 1 - y_eps]` (in `float`), and returns fallback on nonfinite `r`/`t` or failed transforms/prediction.
- **Primitive completion (`c2p_nn_guess_primitives.c`):** requires a validated model, positive finite `D`, finite auxiliaries, and ordered bounds. On valid candidates, clamps `x` to `[1 + q - s, 2 + 2q - s]` and calls `ghl_tabulated_primitive_guess_from_x`. On missing models or unusable states, returns atmosphere fields (`vU = -betaU`, `u0 = lapseinv`, preserved `BU`, and finite conservative `Y_e` clamped to EOS bounds); if only EOS inversion fails after a valid algebraic candidate, keeps bounded `(rho, Y_e, eps)`, `T_max`, and atmosphere `press`/`entropy` while completing velocity.

Evidence (descriptive behavior): `c2p_nn_guess_x.c::ghl_c2p_nn_guess`; corroboration: `ghl_con2prim.h::GHL_NN_C2P_API_VERSION` and `c2p_nn_guess_primitives.c`.

### HDF5 Model Layout, Validation, Build Rules, And Tests

- **HDF5 locations (`c2p_nn_load_from_eos_hdf5.c`):** `ghl_c2p_nn_load_hdf5` reads from file root; `ghl_c2p_nn_load_from_eos_hdf5` reads from group `grhayl_nn_c2p`. Datasets: `dims/{in_dim,hidden_dim,n_hidden,out_dim}` (each `1..8`), `meta/{q_idx,s_idx,y_eps,dx_eps}`, `scaling/{x_eps,x_kind,x_lo,x_hi,x_invrng,out_kind,out_lo,out_hi,out_invrng}`, and `layers/{W_in,b_in,W_hid,b_hid,W_out,b_out}` (`W_hid`/`b_hid` required only for `n_hidden > 1`). Legacy single-output files without `scaling/out_kind` (`out_dim == 1`) default to bounded `x` on `[0, 1]`.
- **Validation (`c2p_nn_validate_model.c`):** enforces `in_dim == 4`, `q_idx == 0`, `s_idx == 2` (incompatible older models must be re-exported, e.g. CI uses `pyghl==0.4.0`), positive dims, finite positive `x_eps`/`dx_eps`, `0 < y_eps < 0.5`, input kinds in `{0, 1}`, output `0` kind `GHL_NN_C2P_OUT_X_BOUNDED`, positive ranges, and finite weights/biases. Successful loads replace and free old `eos->c2p_nn`; failed loads preserve existing `eos->c2p_nn`. `NRPyEOS_free_memory.c` calls `ghl_c2p_nn_free` (requires heap-allocated model/arrays) and clears `eos->c2p_nn = NULL`.
- **Build, retry, GRHayLib & tests:** `neural_network_guess/make.code.defn` registers `#! INCS = ghl_c2p_nn.h` for dependency rebuilds (not installed). Under `GHL_DISABLE_HDF5`, inference/validation remain compiled while HDF5 loaders return `ghl_error_used_disabled_hdf5`. `GRHayLib` maps `enable_backup_nn_primitive_guess` (default `"no"`) to `enable_neural_net_c2p`, treating incompatible embedded metadata as fatal at init. Tested by `unit_test_c2p_nn_guess.c`, `unit_test_con2prim_tabulated.c` (with `pyghl append`), and `unit_test_code_error.c` keys `83`–`85`.

Evidence (descriptive behavior): `c2p_nn_validate_model.c::ghl_c2p_nn_validate_model`; corroboration: `c2p_nn_load_from_eos_hdf5.c` and `NRPyEOS_free_memory.c`.

## Sources

- [GRHayL/include/ghl_con2prim.h](../../../GRHayL/include/ghl_con2prim.h) - `GHL_NN_C2P_API_VERSION`, `ghl_c2p_nn_model`, NN API
- [GRHayL/include/ghl.h](../../../GRHayL/include/ghl.h) - `enable_neural_net_c2p`, `c2p_nn`, `ghl_error_nn_c2p_*`
- [GRHayL/Con2Prim/Tabulated/tabulated_primitive_guess_helpers.c](../../../GRHayL/Con2Prim/Tabulated/tabulated_primitive_guess_helpers.c) - `{q, r, s, t}` auxiliaries and completion helper
- [GRHayL/Con2Prim/Tabulated/neural_network_guess/ghl_c2p_nn.h](../../../GRHayL/Con2Prim/Tabulated/neural_network_guess/ghl_c2p_nn.h) - private inline activation/scaling helpers
- [GRHayL/Con2Prim/Tabulated/neural_network_guess/c2p_nn_guess_x.c](../../../GRHayL/Con2Prim/Tabulated/neural_network_guess/c2p_nn_guess_x.c) - `ghl_c2p_nn_guess` and midpoint fallback
- [GRHayL/Con2Prim/Tabulated/neural_network_guess/c2p_nn_guess_primitives.c](../../../GRHayL/Con2Prim/Tabulated/neural_network_guess/c2p_nn_guess_primitives.c) - `ghl_c2p_nn_guess_primitives` and atmosphere fallback
- [GRHayL/Con2Prim/Tabulated/neural_network_guess/c2p_nn_load_from_eos_hdf5.c](../../../GRHayL/Con2Prim/Tabulated/neural_network_guess/c2p_nn_load_from_eos_hdf5.c) - root/embedded HDF5 loaders and no-HDF5 stubs
- [GRHayL/Con2Prim/Tabulated/neural_network_guess/c2p_nn_validate_model.c](../../../GRHayL/Con2Prim/Tabulated/neural_network_guess/c2p_nn_validate_model.c) - dimension and metadata checks
- [GRHayL/Con2Prim/Tabulated/neural_network_guess/c2p_nn_free.c](../../../GRHayL/Con2Prim/Tabulated/neural_network_guess/c2p_nn_free.c) - heap cleanup
- [GRHayL/Con2Prim/Tabulated/neural_network_guess/make.code.defn](../../../GRHayL/Con2Prim/Tabulated/neural_network_guess/make.code.defn) - `#! INCS = ghl_c2p_nn.h`
- [GRHayL/EOS/Tabulated/NRPyEOS_free_memory.c](../../../GRHayL/EOS/Tabulated/NRPyEOS_free_memory.c) - releases `eos->c2p_nn`
- [Unit_Tests/unit_test_c2p_nn_guess.c](../../../Unit_Tests/unit_test_c2p_nn_guess.c) - direct NN validation/inference/HDF5 tests

## See Also

- Depends on: [Recovery Flow](recovery-flow.md)
- Depends on: [Atmosphere](../atmosphere.md)
- See also: [GRHayLib Runtime Parameters](../../implementations/grhaylib-runtime-parameters.md)
- Validated by: [Con2Prim Tests And Reference Data](tests-and-reference-data.md)
