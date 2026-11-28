# Errors, IO, Debug, And Utilities

> Core `ghl_abort_if_error` switch table, `ghl_Warn_Error` termination rules, `ghl_debug.h` print helpers, and `min_max_clamp.c` utilities. · Status: confirmed
> Up: [GRHayL Core](index.md)

## Summary

Core separates recoverable `ghl_error_codes_t` return values from `ghl_io.h` reporting macros (`ghl_warn` returns; `ghl_error`, `ghl_abort`, and nonzero-code `ghl_Error` terminate), maps error codes to messages in `abort_if_error.c`, supplies stderr debug printers in `ghl_debug.h`, and provides unchecked integer/floating-point clamp helpers in `min_max_clamp.c`.

## Detail

### Error Codes, `ghl_abort_if_error`, And Terminating IO

- **`ghl_abort_if_error` (`abort_if_error.c`):** switches over every declared `ghl_error_codes_t` value in `ghl.h`: `ghl_success` returns normally, while known non-success codes call `ghl_Error(error, ...)` which invokes `exit(error)` (not `abort()`) because current enum codes are positive. Because there is no `default` switch arm, an out-of-range cast enum returns silently without terminating.
- **`ghl_info` and `ghl_Warn_Error` (`info_warn_error.c`, `ghl_io.h`):** `ghl_info` prints a GRHayL-prefixed message to stdout. `ghl_warn`, `ghl_error`, `ghl_abort`, and `ghl_Error` pass `__FILE__`, `__LINE__`, and `__func__` to `ghl_Warn_Error`, which terminates strictly by `exit_code`: `0` returns (`ghl_warn`), `-1` (`ghl_error_abort`) calls `abort()` (`ghl_abort`), and any other nonzero value calls `exit(exit_code)` (`ghl_error` exits `1`; `ghl_Error` uses the caller's code).
- **Expected-failure testing:** [Expected Failure And Error Keys](../tests/expected-failure-and-error-keys.md) describes how `unit_test_code_error.c` keys exit: matched error codes through `ghl_abort_if_error`, and other expected outcomes, including no-HDF5 skips, through `pass_test`.

Evidence (descriptive behavior): `abort_if_error.c::ghl_abort_if_error`; corroboration: `info_warn_error.c::ghl_Warn_Error`.

### Debug Printers And Min/Max/Clamp Utilities

- **`ghl_debug.h`:** defines inline stderr printers with function/line wrapper macros. `ghl_debug_print_prims` prints thermodynamic fields, `vU[0..2]`, and `BU[0..2]` while omitting `u0`; `ghl_debug_print_cons` prints `D`/`rho`, `tau`, `SD0`, `SD1`, `SD2`, `entropy`, and `Y_e` in that order. `unit_test_con2prim_debug.c` invokes both without asserting printed text.
- **`min_max_clamp.c`:** implements `ghl_imin`, `ghl_imax`, `ghl_iclamp`, and `ghl_clamp`. Tabulated EOS interpolators instead call a private always-inline `NRPyEOS_iclamp` in `NRPyEOS_tabulated_helpers.h` with identical `min(max(x, x_min), x_max)` order—keep both synchronized if clamp semantics change. Neither clamp validates `x_min <= x_max`, and no direct unit test asserts Core clamp functions.

Evidence (descriptive behavior): `min_max_clamp.c::ghl_iclamp`; corroboration: `NRPyEOS_tabulated_helpers.h::NRPyEOS_iclamp` and `ghl_debug.h::ghl_debug_print_prims`.

## Sources

- [GRHayL/include/ghl.h](../../GRHayL/include/ghl.h) - `ghl_error_codes_t`, `ghl_abort_if_error`, clamp declarations
- [GRHayL/include/ghl_io.h](../../GRHayL/include/ghl_io.h) - `ghl_info`, `ghl_warn`, `ghl_error`, `ghl_abort`, `ghl_Error`
- [GRHayL/include/ghl_debug.h](../../GRHayL/include/ghl_debug.h) - `ghl_debug_print_prims`, `ghl_debug_print_cons`
- [GRHayL/GRHayL_Core/abort_if_error.c](../../GRHayL/GRHayL_Core/abort_if_error.c) - error-code message switch
- [GRHayL/GRHayL_Core/info_warn_error.c](../../GRHayL/GRHayL_Core/info_warn_error.c) - `ghl_info` and `ghl_Warn_Error`
- [GRHayL/GRHayL_Core/min_max_clamp.c](../../GRHayL/GRHayL_Core/min_max_clamp.c) - `ghl_imin`, `ghl_imax`, `ghl_iclamp`, `ghl_clamp`
- [GRHayL/EOS/Tabulated/interpolators/NRPyEOS_tabulated_helpers.h](../../GRHayL/EOS/Tabulated/interpolators/NRPyEOS_tabulated_helpers.h) - private inline `NRPyEOS_iclamp`

## See Also

- Depends on: [Shared Parameters And Enums](shared-parameters-and-enums.md)
- Validated by: [Expected Failure And Error Keys](../tests/expected-failure-and-error-keys.md)
- See also: [Core Tests And Reference Data](tests-and-reference-data.md)
