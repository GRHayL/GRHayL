# Expected Failure And Error Keys

> Inverted exit-code protocol, `GHL_DISABLE_HDF5` skip sets, and subsystem key ranges (`0`–`90`) in `Unit_Tests/unit_test_code_error.c`. · Status: confirmed
> Up: [Tests](index.md)

## Summary

`Unit_Tests/unit_test_code_error.c` verifies GRHayL error-code paths over integer keys `0`–`90` using an inverted process-status protocol: matched errors exit nonzero through `ghl_abort_if_error`, whereas unexpected internal setup failures exit `0` so `.github/run_tests.sh` reports `Failed to fail!`. Under `GHL_DISABLE_HDF5`, HDF5-dependent keys exit as expected skips rather than running tabulated code.

## Detail

### Inverted Exit Protocol And `GHL_DISABLE_HDF5` Skips

Claim status: contested; contradiction: CONTR-0010. See [CONTR-0010](../contradictions.md#contr-0010).

`.github/run_tests.sh` runs `unit_test_code_error` over keys `0..90` inside a private `code_error_workdir` removed on `EXIT`. Inside `unit_test_code_error.c`:
- `expect_error_code(...)` compares the produced `ghl_error_codes_t` against the expected code, then invokes `ghl_abort_if_error(error)` for nonzero termination (the test does not assert emitted message strings).
- `pass_test(...)` exits nonzero for intended expected-error or skip outcomes; `fail_test(...)` exits `0` on unexpected internal setup failures (such as common hybrid setup failure, key `33` returning non-`NULL`, or key `80` prerequisite failure) so the outer runner fails the key.
- Keys `86`–`88` pass enum values directly (`ghl_abort_if_error` coverage for NRPyLeakage blocking/nonfinite errors and `ghl_error_invalid_hlle_wavespeeds`) rather than invoking the producing physics routines.
- Under `GHL_DISABLE_HDF5`, keys `6..32`, `34..60`, `63`, `66`, `69..85`, and `89..90` take expected-skip exits.

Evidence (descriptive behavior): `Unit_Tests/unit_test_code_error.c::expect_error_code`; corroboration: `Unit_Tests/unit_test_code_error.c::fail_test` and `.github/run_tests.sh::Failed to fail!`.

### Error-Key Ranges By Subsystem

| Key range | Subsystem & behavior checked | Defining route |
| --- | --- | --- |
| `0`, `1`, `61`, `62`, `64`, `65`, `67`, `68` | Invalid simple/hybrid EOS density and pressure setup. | [EOS Dispatch](../core/eos-dispatch.md) |
| `2`, `3` | Invalid `NRPyLeakage_Fermi_Dirac_integrals` keys across the small-`z` split. | [Neutrinos Tests](../gems/neutrinos/tests-and-reference-data.md) |
| `4`, `5`, `32`, `33` | Core `ghl_limit_v_and_compute_u0` singular `u0` (`4`); invalid Con2Prim routine selection (`5`, `32`) and `NULL` routine-name lookup (`33`). | [Velocity And `u0`](../core/velocity-u0.md), [Recovery Flow](../gems/con2prim/recovery-flow.md) |
| `6..11`, `63`, `66`, `69..72` | Invalid tabulated `rho`, `Y_e`, `T` atmosphere/min/max bounds (`SLy4` table). | [Tabulated Table Layout](../gems/eos/tabulated-table-layout.md) |
| `12..31` | Too-many-variable helper calls and out-of-table `rho`/`Y_e`/`T` bounds in direct and root-finding interpolators. | [Tabulated Interpolation And Bounds](../gems/eos/tabulated-interpolation-and-bounds.md) |
| `34..60`, `73`, `74` | Missing (`34`) and malformed temporary `test.h5` HDF5 open/rank/read failures (standalone invocations can leave `test.h5` in cwd). | [HDF5 Sample Tables](hdf5-sample-tables.md) |
| `75..77` | Null `ghl_eos_parameters`, non-tabulated EOS type, and invalid `ghl_eos_table_t`. | [EOS Initialization And Dispatch](../gems/eos/initialization-and-dispatch.md) |
| `78..82` | Beta-equilibrium out-of-range density/pressure `root not bracketed` errors. | [Tabulated Interpolation And Bounds](../gems/eos/tabulated-interpolation-and-bounds.md) |
| `83..85` | `ghl_c2p_nn_load_hdf5` null EOS (`83`), empty path (`84`), and `enable_neural_net_c2p` on `SLy4` lacking `grhayl_nn_c2p` (`85`). | [Neural-Network Primitive Guess](../gems/con2prim/neural-network-primitive-guess.md) |
| `86..88` | Direct enum abort checks for NRPyLeakage blocking/nonfinite errors (`86`, `87`) and `ghl_error_invalid_hlle_wavespeeds` (`88`). | [Errors, IO, Debug, And Utilities](../core/errors-io-debug-utilities.md) |
| `89`, `90` | Tabulated `ghl_tabulated_Newman1D_energy` (`89`) and `ghl_tabulated_Palenzuela1D_energy` (`90`) returning `ghl_error_table_bisection` on over-range `tau`. | [Recovery Flow](../gems/con2prim/recovery-flow.md) |

## Sources

- [Unit_Tests/unit_test_code_error.c](../../Unit_Tests/unit_test_code_error.c) - keys `0`–`90`, `expect_error_code`, `pass_test`, `fail_test`, `GHL_DISABLE_HDF5` skips
- [.github/run_tests.sh](../../.github/run_tests.sh) - `code_error_workdir` loop and `Failed to fail!` check

## See Also

- Depends on: [Errors, IO, Debug, And Utilities](../core/errors-io-debug-utilities.md)
- Depends on: [Runner And Build Outputs](runner-and-build-outputs.md)
- See also: [HDF5 Sample Tables](hdf5-sample-tables.md)
- See also: [Contradictions](../contradictions.md#contr-0010)
