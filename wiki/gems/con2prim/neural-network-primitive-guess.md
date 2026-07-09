# Neural-Network Primitive Guess

This page routes the tabulated Con2Prim neural-network primitive guess feature. It covers the GRHayL API, embedded HDF5 model contract, runtime retry flow, tests, CI setup, and downstream GRHayLib flag. It does not document a standalone solver, a `ghl_con2prim_id_t` value, hybrid/simple EOS support, model training, or model performance claims.

## Purpose

The neural-network feature predicts an initial scalar `x` used to rebuild a tabulated primitive guess. Normal tabulated Con2Prim routines still perform the recovery; the NN path only changes the guess used for retry after a failure. The owning runtime path is tabulated Con2Prim in [GRHayL/Con2Prim/con2prim_multi_method.c](../../../GRHayL/Con2Prim/con2prim_multi_method.c), with helper completion in [GRHayL/Con2Prim/Tabulated/tabulated_primitive_guess_helpers.c](../../../GRHayL/Con2Prim/Tabulated/tabulated_primitive_guess_helpers.c).

It is not a solver enum: [GRHayL/include/ghl.h](../../../GRHayL/include/ghl.h) lists `ghl_con2prim_id_t` routines, while the NN structs and functions are a separate "Neural Network Initial Guesses" API in [GRHayL/include/ghl_con2prim.h](../../../GRHayL/include/ghl_con2prim.h).

## Source Authority

Primary source paths:

- Public API: [GRHayL/include/ghl_con2prim.h](../../../GRHayL/include/ghl_con2prim.h) and [GRHayL/include/ghl.h](../../../GRHayL/include/ghl.h).
- NN implementation: [GRHayL/Con2Prim/Tabulated/neural_network_guess/](../../../GRHayL/Con2Prim/Tabulated/neural_network_guess/).
- EOS initialization: [GRHayL/GRHayL_Core/initialize_eos.c](../../../GRHayL/GRHayL_Core/initialize_eos.c).
- Runtime retry flow: [GRHayL/Con2Prim/con2prim_multi_method.c](../../../GRHayL/Con2Prim/con2prim_multi_method.c).
- Primitive completion: [GRHayL/Con2Prim/Tabulated/tabulated_primitive_guess_helpers.c](../../../GRHayL/Con2Prim/Tabulated/tabulated_primitive_guess_helpers.c).
- Lifetime/free path: [GRHayL/EOS/Tabulated/NRPyEOS_free_memory.c](../../../GRHayL/EOS/Tabulated/NRPyEOS_free_memory.c).
- Build filtering: [configure](../../../configure), [GRHayL/Con2Prim/Tabulated/make.code.defn](../../../GRHayL/Con2Prim/Tabulated/make.code.defn), and [GRHayL/Con2Prim/Tabulated/neural_network_guess/make.code.defn](../../../GRHayL/Con2Prim/Tabulated/neural_network_guess/make.code.defn).
- Tests and CI: [Unit_Tests/unit_test_c2p_nn_guess.c](../../../Unit_Tests/unit_test_c2p_nn_guess.c), [Unit_Tests/unit_test_con2prim_tabulated.c](../../../Unit_Tests/unit_test_con2prim_tabulated.c), [Unit_Tests/unit_test_code_error.c](../../../Unit_Tests/unit_test_code_error.c), [.github/run_tests.sh](../../../.github/run_tests.sh), and [.github/workflows/](../../../.github/workflows/).
- Downstream flag routing: [implementations/GRHayLib/param.ccl](../../../implementations/GRHayLib/param.ccl) and [implementations/GRHayLib/src/initialize_and_shutdown.c](../../../implementations/GRHayLib/src/initialize_and_shutdown.c).

## Public API

[GRHayL/include/ghl_con2prim.h](../../../GRHayL/include/ghl_con2prim.h) defines:

- `ghl_nn_c2p_input_t`: four `float` inputs `q`, `r`, `s`, and `t`.
- `ghl_nn_c2p_guess_t`: one `float x` prediction.
- `ghl_c2p_nn_model`: dimensions, input indices, numeric epsilons, scaling arrays, output metadata, and layer weights/biases.
- `GHL_NN_C2P_API_VERSION`: public schema/output-kind version, currently `3u`.
- `ghl_c2p_nn_guess`: maps `ghl_nn_c2p_input_t` plus `ghl_c2p_nn_model` to an `x` guess.
- `ghl_c2p_nn_guess_primitives`: computes tabulated auxiliaries, asks the NN for `x`, clamps it, and completes primitives.
- `ghl_c2p_nn_validate_model`: validates dimensions, indices, scaling metadata, output kinds, required arrays, and finite weights.
- `ghl_c2p_nn_free`: releases all model arrays and the model struct.
- `ghl_c2p_nn_load_hdf5`: loads a standalone/root HDF5 NN model into `eos->c2p_nn`.
- `ghl_c2p_nn_load_from_eos_hdf5`: loads the embedded EOS-table group into `eos->c2p_nn`.

[GRHayL/include/ghl.h](../../../GRHayL/include/ghl.h) carries the shared state and errors: `ghl_eos_parameters.enable_neural_net_c2p`, `ghl_eos_parameters.c2p_nn`, and NN-specific errors such as `ghl_error_nn_c2p_model_is_null`, `ghl_error_nn_c2p_invalid_dimensions`, `ghl_error_nn_c2p_invalid_input_index`, `ghl_error_nn_c2p_missing_array`, `ghl_error_nn_c2p_invalid_kind`, and `ghl_error_nn_c2p_invalid_number`.

## Runtime Flow

[GRHayL/GRHayL_Core/initialize_eos.c](../../../GRHayL/GRHayL_Core/initialize_eos.c) sets `eos->enable_neural_net_c2p` from `ghl_initialize_tabulated_eos(... enable_neural_net_c2p ...)`, initializes `eos->c2p_nn` to `NULL`, reads the tabulated EOS table, and, when the NN flag is true, calls `ghl_c2p_nn_load_from_eos_hdf5(table_filepath, eos)`. If that load fails, initialization frees the partially initialized tabulated EOS memory and returns the error.

[GRHayL/Con2Prim/con2prim_multi_method.c](../../../GRHayL/Con2Prim/con2prim_multi_method.c) runs the usual `params->calc_prim_guess` path first, stores that primitive guess, and calls the selected tabulated main routine. If the main routine fails and `eos->enable_neural_net_c2p` is true, it builds `prims_guess_nn` from the stored guess, calls `ghl_c2p_nn_guess_primitives`, and retries the same main routine with the NN-derived guess.

The backup loop in [GRHayL/Con2Prim/con2prim_multi_method.c](../../../GRHayL/Con2Prim/con2prim_multi_method.c) first resets primitives to the ordinary stored guess and runs each configured tabulated backup. If that backup fails and the NN flag is true, it retries the same backup with `prims_guess_nn`. This makes the NN path a failure retry policy for tabulated Con2Prim only.

## Numerical Contract

[GRHayL/Con2Prim/Tabulated/tabulated_primitive_guess_helpers.c](../../../GRHayL/Con2Prim/Tabulated/tabulated_primitive_guess_helpers.c) computes shared auxiliaries from conservative variables and the metric: `q = tau / D`, `r = S_squared / D^2`, `s = B_squared / D`, and `t = BdotS / D^(3/2)`, where `D` is `cons_undens->rho`.

[GRHayL/Con2Prim/Tabulated/neural_network_guess/c2p_nn_guess_x.c](../../../GRHayL/Con2Prim/Tabulated/neural_network_guess/c2p_nn_guess_x.c) uses `q`, `r`, `s`, and `t` as NN inputs. If the model is `NULL`, `q` or `s` is non-finite, or the derived `x_lo = 1 + q - s` or width `1 + q` is non-finite, it returns the zero-initialized guess. Otherwise it sets fallback `x` to the midpoint of `x_lo` and the width, with width floored by `model->dx_eps`; that midpoint fallback is returned when `r` or `t` is non-finite, feature transforms fail, prediction fails, or the predicted `x` is non-finite.

[GRHayL/Con2Prim/Tabulated/neural_network_guess/ghl_c2p_nn.h](../../../GRHayL/Con2Prim/Tabulated/neural_network_guess/ghl_c2p_nn.h) transforms input features by kind `0` as identity and kind `1` as `log10` after an `x_eps` floor, clips transformed inputs into the configured scaling interval, applies hard-tanh hidden layers, applies sigmoid output activation, and clamps each output into `(y_eps, 1 - y_eps)`.

[GRHayL/Con2Prim/Tabulated/neural_network_guess/c2p_nn_guess_primitives.c](../../../GRHayL/Con2Prim/Tabulated/neural_network_guess/c2p_nn_guess_primitives.c) clamps the guessed `x` to `[1 + q - s, 2 + 2q - s]`, then calls `ghl_tabulated_primitive_guess_from_x`. That helper in [GRHayL/Con2Prim/Tabulated/tabulated_primitive_guess_helpers.c](../../../GRHayL/Con2Prim/Tabulated/tabulated_primitive_guess_helpers.c) computes Lorentz factor, `rho`, `Y_e`, `eps`, `u0`, pressure, entropy, temperature, and velocity through existing tabulated EOS and speed-limit helpers.

## HDF5 Model Contract

[GRHayL/Con2Prim/Tabulated/neural_network_guess/c2p_nn_load_from_eos_hdf5.c](../../../GRHayL/Con2Prim/Tabulated/neural_network_guess/c2p_nn_load_from_eos_hdf5.c) supports two load locations:

- Root model files through `ghl_c2p_nn_load_hdf5`, where datasets are read directly from the file root.
- Embedded EOS-table models through `ghl_c2p_nn_load_from_eos_hdf5`, where the group prefix must be `grhayl_nn_c2p`.

The loader reads scalar datasets under `dims/` for `in_dim`, `hidden_dim`, `n_hidden`, and `out_dim`; scalar datasets under `meta/` for `q_idx`, `s_idx`, `y_eps`, and `dx_eps`; scalar dataset `scaling/x_eps`; array datasets under `scaling/` for `x_kind`, `x_lo`, `x_hi`, `x_invrng`, `out_kind`, `out_lo`, `out_hi`, and `out_invrng`; and layer arrays under `layers/` for `W_in`, `b_in`, optional `W_hid`, optional `b_hid`, `W_out`, and `b_out`.

The loader allows dimensions from `1` through `8` before allocating arrays, then calls `ghl_c2p_nn_validate_model`. Validation in [GRHayL/Con2Prim/Tabulated/neural_network_guess/c2p_nn_validate_model.c](../../../GRHayL/Con2Prim/Tabulated/neural_network_guess/c2p_nn_validate_model.c) requires `in_dim == 4`, positive hidden/output dimensions, valid `q_idx` and `s_idx`, positive numeric epsilons, required arrays, input kinds `0` or `1`, output kind `GHL_NN_C2P_OUT_X_BOUNDED` for output `0`, linear or log-linear kinds for later outputs, positive ranges, and finite weights/biases.

Legacy one-output files are accepted when `scaling/out_kind` is absent and `out_dim == 1`; the loader fills output scaling as bounded `x` on `[0, 1]`. Missing output scaling with more than one output is an HDF5 dataset-open error. These paths are exercised in [Unit_Tests/unit_test_c2p_nn_guess.c](../../../Unit_Tests/unit_test_c2p_nn_guess.c).

Successful loads free any old `eos->c2p_nn` and replace it with the new model. Failed loads free the new partial model and preserve any existing `eos->c2p_nn`, as implemented in [GRHayL/Con2Prim/Tabulated/neural_network_guess/c2p_nn_load_from_eos_hdf5.c](../../../GRHayL/Con2Prim/Tabulated/neural_network_guess/c2p_nn_load_from_eos_hdf5.c) and tested by `check_hdf5_load_failure_preserves_model` in [Unit_Tests/unit_test_c2p_nn_guess.c](../../../Unit_Tests/unit_test_c2p_nn_guess.c).

## Build/Lifetime

[GRHayL/Con2Prim/Tabulated/make.code.defn](../../../GRHayL/Con2Prim/Tabulated/make.code.defn) includes `neural_network_guess`, and [GRHayL/Con2Prim/Tabulated/neural_network_guess/make.code.defn](../../../GRHayL/Con2Prim/Tabulated/neural_network_guess/make.code.defn) lists the NN free, guess, loader, and validation sources.

When HDF5 is disabled, [GRHayL/Con2Prim/Tabulated/neural_network_guess/c2p_nn_load_from_eos_hdf5.c](../../../GRHayL/Con2Prim/Tabulated/neural_network_guess/c2p_nn_load_from_eos_hdf5.c) compiles loader stubs returning `ghl_error_used_disabled_hdf5`. [configure](../../../configure) defines `GHL_DISABLE_HDF5` and filters tabulated source lists while retaining the tabulated primitive-guess helper and `neural_network_guess` sources.

[GRHayL/EOS/Tabulated/NRPyEOS_free_memory.c](../../../GRHayL/EOS/Tabulated/NRPyEOS_free_memory.c) frees `eos->c2p_nn` with `ghl_c2p_nn_free` and then sets it to `NULL` when tabulated EOS memory is released.

## Tests/CI

[Unit_Tests/unit_test_c2p_nn_guess.c](../../../Unit_Tests/unit_test_c2p_nn_guess.c) covers model validation errors, direct `ghl_c2p_nn_guess`, fallback behavior for invalid/non-finite inputs, root HDF5 loading, embedded `grhayl_nn_c2p` loading, legacy one-output loading, malformed HDF5 datasets, validation failure after load, and preservation of an existing model after a failed load.

[Unit_Tests/unit_test_con2prim_tabulated.c](../../../Unit_Tests/unit_test_con2prim_tabulated.c) enables embedded NN loading when the run key is nonzero, then replays tabulated Con2Prim tests with `eos.enable_neural_net_c2p` both disabled and enabled. Its generation path explicitly disables NN guesses.

[Unit_Tests/unit_test_code_error.c](../../../Unit_Tests/unit_test_code_error.c) covers NN-related expected-failure keys `83`, `84`, and `85`: null EOS pointer for `ghl_c2p_nn_load_hdf5`, empty NN filepath for `ghl_c2p_nn_load_hdf5`, and NN-enabled tabulated EOS initialization against a table without embedded `grhayl_nn_c2p`.

[.github/run_tests.sh](../../../.github/run_tests.sh) runs `./test/unit_test_c2p_nn_guess`, loops `unit_test_code_error` over `0..85`, and runs `pyghl append SLy4_3335_rho391_temp163_ye66.h5` before tabulated Con2Prim replay. The workflow files in [.github/workflows/](../../../.github/workflows/) include `c2p_nn_guess` in the failure-test matrix, run the same `pyghl append` command, and loop error keys over `0..85`. This page treats `pyghl append` only as visible test setup evidence; it does not document `pyghl` internals.

## Downstream Flag

[implementations/GRHayLib/param.ccl](../../../implementations/GRHayLib/param.ccl) defines `enable_backup_nn_primitive_guess`, defaulting to `"no"`, with text saying `"yes"` uses the neural-network primitive guess if the default guess fails. [implementations/GRHayLib/src/initialize_and_shutdown.c](../../../implementations/GRHayLib/src/initialize_and_shutdown.c) passes that Cactus parameter as the `enable_neural_net_c2p` argument to `ghl_initialize_tabulated_eos`.

This downstream flag maps into the tabulated EOS initialization path only. It does not add hybrid/simple EOS support and does not create a new Con2Prim solver ID.

## Known Gaps

No repo-local model training, provenance, or accuracy source was found in the source, tests, CI, Doxygen raw docs, or GRHayLib implementation paths listed above. Do not add claims about those topics unless maintainers add repo-local or official evidence.
