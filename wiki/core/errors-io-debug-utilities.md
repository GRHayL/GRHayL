# Core Errors, IO, Debug, And Utilities

## Purpose

This page routes Core error-code handling, IO helpers, debug print helpers, and
small utility functions. It is not a replacement for headers or source; use it
to find the owning files and tests before changing behavior.

Ground truth:
[`GRHayL/include/ghl.h`](../../GRHayL/include/ghl.h),
[`GRHayL/include/ghl_io.h`](../../GRHayL/include/ghl_io.h),
[`GRHayL/include/ghl_debug.h`](../../GRHayL/include/ghl_debug.h),
[`GRHayL/GRHayL_Core/abort_if_error.c`](../../GRHayL/GRHayL_Core/abort_if_error.c),
[`GRHayL/GRHayL_Core/info_warn_error.c`](../../GRHayL/GRHayL_Core/info_warn_error.c),
and
[`GRHayL/GRHayL_Core/min_max_clamp.c`](../../GRHayL/GRHayL_Core/min_max_clamp.c).

## Error Codes And Message Mapping

`ghl_error_codes_t` lives in
[`GRHayL/include/ghl.h`](../../GRHayL/include/ghl.h). Use that enum as the
authority for available return codes, including EOS, Con2Prim, table, HDF5,
atmosphere, and utility error keys.

`ghl_abort_if_error` is declared in
[`GRHayL/include/ghl.h`](../../GRHayL/include/ghl.h) and implemented in
[`GRHayL/GRHayL_Core/abort_if_error.c`](../../GRHayL/GRHayL_Core/abort_if_error.c).
That source owns the mapping from `ghl_error_codes_t` values to human-readable
messages; do not copy the full table into KB pages.

Every currently declared `ghl_error_codes_t` value has a switch case, and
`ghl_success` returns normally. Known non-success values route through
`ghl_Error(error, ...)`, which calls `exit(error)` because current error codes
are positive. Despite the helper name, these cases do not use `abort()`.
There is no `default` switch arm: a cast/out-of-range value not represented by
the enum returns silently without terminating
([`GRHayL/GRHayL_Core/abort_if_error.c`](../../GRHayL/GRHayL_Core/abort_if_error.c),
[`GRHayL/GRHayL_Core/info_warn_error.c`](../../GRHayL/GRHayL_Core/info_warn_error.c)).

`Unit_Tests/unit_test_code_error.c` exercises expected return-code paths and
then calls `ghl_abort_if_error(error)` to verify the process-level failure path
for configured error cases. Source:
[`Unit_Tests/unit_test_code_error.c`](../../Unit_Tests/unit_test_code_error.c).
The test deliberately exits nonzero after a matched error reaches
`ghl_abort_if_error`; [`.github/run_tests.sh`](../../.github/run_tests.sh)
loops keys 0 through 85 and treats a zero process status as "Failed to fail."
This is expected-failure process coverage, not an ordinary success-status test.

## Return Codes Versus Terminating IO

Most recoverable routines return `ghl_error_codes_t` to callers. The IO macros
in [`GRHayL/include/ghl_io.h`](../../GRHayL/include/ghl_io.h) are different:
`ghl_warn` returns after printing, while `ghl_error`, `ghl_abort`, and
`ghl_Error` terminate through `exit` or `abort` according to the exit code
handled in
[`GRHayL/GRHayL_Core/info_warn_error.c`](../../GRHayL/GRHayL_Core/info_warn_error.c).

`ghl_info` writes a GRHayL-prefixed message to stdout. Source:
[`GRHayL/GRHayL_Core/info_warn_error.c`](../../GRHayL/GRHayL_Core/info_warn_error.c).

`ghl_warn`, `ghl_error`, `ghl_abort`, and `ghl_Error` route through
`ghl_Warn_Error` with `__FILE__`, `__LINE__`, and `__func__` context supplied by
the macros in [`GRHayL/include/ghl_io.h`](../../GRHayL/include/ghl_io.h).
`ghl_Warn_Error` prints that file/line/function context and the formatted
message in
[`GRHayL/GRHayL_Core/info_warn_error.c`](../../GRHayL/GRHayL_Core/info_warn_error.c).

Exact termination control is the `exit_code`: zero returns, `-1`
(`ghl_error_abort`) calls `abort()`, and every other value calls
`exit(exit_code)`. Thus `ghl_warn` returns, `ghl_error` exits with 1,
`ghl_abort` aborts, and `ghl_Error` follows its caller-supplied value.

## Debug Helpers

`ghl_debug.h` defines inline helpers for printing primitive and conservative
struct contents to stderr, plus wrapper macros that add the current function and
line. Source: [`GRHayL/include/ghl_debug.h`](../../GRHayL/include/ghl_debug.h).
Installed public header routing is listed in
[`GRHayL/include/make.code.defn`](../../GRHayL/include/make.code.defn).

Current debug output has two bounded hazards:

- `ghl_debug_print_prims` omits `u0`; it prints thermodynamic fields,
  `vU[0..2]`, and `BU[0..2]` only.
- `ghl_debug_print_cons` labels its first five columns `D`, `tau`, `SD0`,
  `SD1`, `SD2`, but passes values in the order `tau`, `SD0`, `SD1`, `SD2`,
  `rho`. Do not use those labels as trustworthy field identification until
  source is corrected.

`Unit_Tests/unit_test_con2prim_debug.c` calls both inline printers for manual
diagnostic output but does not assert column labels or values. No focused debug
format test is visible in the listed Core evidence.

## Min, Max, And Clamp Utilities

Core declares `ghl_imin`, `ghl_imax`, `ghl_iclamp`, and `ghl_clamp` in
[`GRHayL/include/ghl.h`](../../GRHayL/include/ghl.h) and implements them in
[`GRHayL/GRHayL_Core/min_max_clamp.c`](../../GRHayL/GRHayL_Core/min_max_clamp.c).
`ghl_iclamp` composes the integer min/max helpers; `ghl_clamp` uses the math
library floating-point min/max route. Source:
[`GRHayL/GRHayL_Core/min_max_clamp.c`](../../GRHayL/GRHayL_Core/min_max_clamp.c).

Neither clamp validates that `x_min <= x_max`; reversed bounds remain caller
error and produce the literal composed min/max result. These functions expose
no error return.

Direct clamp test coverage is weak from the listed Core task sources: the Core
suite and code-error test do not directly assert `ghl_imin`, `ghl_imax`,
`ghl_iclamp`, or `ghl_clamp` behavior. Broader repo code uses the clamp helpers
in EOS and Con2Prim paths, but that is indirect coverage. Sources:
[`Unit_Tests/unit_test_grhayl_core_test_suite.c`](../../Unit_Tests/unit_test_grhayl_core_test_suite.c),
[`Unit_Tests/unit_test_code_error.c`](../../Unit_Tests/unit_test_code_error.c),
[`GRHayL/GRHayL_Core/min_max_clamp.c`](../../GRHayL/GRHayL_Core/min_max_clamp.c).

## Source Of Truth

- [`wiki/build-and-ci.md`](../build-and-ci.md)
- [`docs/raw/GRHayL_Core.dox`](../../docs/raw/GRHayL_Core.dox)
- [`GRHayL/include/ghl.h`](../../GRHayL/include/ghl.h)
- [`GRHayL/include/ghl_io.h`](../../GRHayL/include/ghl_io.h)
- [`GRHayL/include/ghl_debug.h`](../../GRHayL/include/ghl_debug.h)
- [`GRHayL/include/make.code.defn`](../../GRHayL/include/make.code.defn)
- [`GRHayL/GRHayL_Core/abort_if_error.c`](../../GRHayL/GRHayL_Core/abort_if_error.c)
- [`GRHayL/GRHayL_Core/info_warn_error.c`](../../GRHayL/GRHayL_Core/info_warn_error.c)
- [`GRHayL/GRHayL_Core/min_max_clamp.c`](../../GRHayL/GRHayL_Core/min_max_clamp.c)
- [`Unit_Tests/unit_test_grhayl_core_test_suite.c`](../../Unit_Tests/unit_test_grhayl_core_test_suite.c)
- [`Unit_Tests/unit_test_code_error.c`](../../Unit_Tests/unit_test_code_error.c)
