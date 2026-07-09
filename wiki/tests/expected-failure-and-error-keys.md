# Expected Failure And Error Keys

Purpose: route `unit_test_code_error.c` expected-failure keys, no-HDF5 skip
semantics, and cross-gem ownership of error-path coverage. This page does not
re-document the error-code API; use the linked Core, EOS, Con2Prim, and
Neutrinos pages for behavior contracts.

Read with [Core tests and fixtures](../core/tests-and-fixtures.md),
[Core errors, IO, debug, utilities](../core/errors-io-debug-utilities.md),
[EOS tests and fixtures](../gems/eos/tests-and-fixtures.md),
[EOS initialization and dispatch](../gems/eos/initialization-and-dispatch.md),
[Con2Prim tests and fixtures](../gems/con2prim/tests-and-fixtures.md),
[Con2Prim recovery flow](../gems/con2prim/recovery-flow.md), and
[Neutrinos tests and fixtures](../gems/neutrinos/tests-and-fixtures.md).

## Harness Contract

The full runner invokes `./test/unit_test_code_error "$i"` for keys `0..85`.
Each key is treated as an expected-error case: if the executable exits
successfully, `.github/run_tests.sh` prints `Failed to fail!` and fails the
runner; if the executable exits nonzero, the runner treats that as the expected
failure and continues.

Inside `unit_test_code_error.c`, `expect_error_code(...)` checks that a call
returned the expected GRHayL error code, then routes through
`ghl_abort_if_error(error)` so the process exits as the expected-error harness
requires. Key `33` is a special invalid-name route: it expects
`ghl_get_con2prim_routine_name` to return `NULL` and exits nonzero through the
same expected-failure convention.

## No-HDF5 Skip Semantics

When compiled with `GHL_DISABLE_HDF5`, `unit_test_code_error.c` skips exactly
these HDF5-only keys:

- `6..32`
- `34..60`
- `63`
- `66`
- `69..85`

Those skipped HDF5-only keys exit as the test harness expects under a no-HDF5
build. Treat them as expected-skip confirmations for HDF5-only paths, not as
ordinary behavior coverage. The build-side no-HDF5 source and test filtering is
owned by [configure](../../configure), which defines `GHL_DISABLE_HDF5` and
filters tabulated/HDF5 tests and sources.

## Behavior Families

- Hybrid/simple EOS parameter validation: keys `0`, `1`, `61`, `62`, `64`,
  `65`, `67`, and `68` cover invalid density and pressure setup for simple and
  hybrid EOS paths. Route contract questions to Core EOS setup and EOS
  initialization pages.
- Neutrinos Fermi-Dirac invalid-key handling: keys `2` and `3` call
  `NRPyLeakage_Fermi_Dirac_integrals` on both sides of the small-`z` split.
  Route behavior to Neutrinos tests and fixtures.
- Core velocity and Con2Prim dispatch: key `4` covers singular `u0` through
  `ghl_limit_v_and_compute_u0`; keys `5`, `32`, and `33` cover invalid
  Con2Prim routine selection or invalid routine-name lookup. Route `u0` to Core
  velocity pages and routine dispatch to Con2Prim recovery flow.
- Tabulated EOS setup bounds: keys `6..11`, `63`, `66`, and `69..72` cover
  invalid tabulated `rho`, `Y_e`, and `T` atmosphere/min/max setup. Route to EOS
  initialization and tabulated table pages.
- Tabulated interpolation and bounds: keys `12..31` cover too-many-variable
  helper calls plus out-of-table `rho`, `Y_e`, and `T` behavior for direct
  `(rho,Y_e,T)` interpolators and root-finding auxiliary routes. These are EOS
  table/interpolator test routes, not a replacement for the public API docs.
- HDF5 table-read failures: key `34` removes any existing `test.h5` and checks
  the missing-file path; keys `35..60`, `73`, and `74` create temporary
  malformed `test.h5` files to cover dataset-open, dataset-rank, and
  dataset-read failures. Route table-shape and HDF5 questions to EOS tests and
  HDF5 sample table routing.
- Invalid tabulated EOS/table state: keys `75..77` cover null EOS parameter
  struct, non-tabulated EOS type, and invalid table type before or during table
  setup. Route to EOS initialization and Core shared parameter pages.
- Beta-equilibrium/root-bracketing failures: keys `78..82` initialize
  beta-equilibrium helpers, then force out-of-range density or pressure routes
  that return root not bracketed. Route detailed root behavior to tabulated EOS
  interpolation/bounds pages.
- NN loader/init failures: key `83` passes a NULL EOS parameter struct to
  `ghl_c2p_nn_load_hdf5`; key `84` passes an empty NN filepath to
  `ghl_c2p_nn_load_hdf5`; key `85` initializes tabulated EOS with
  `enable_neural_net_c2p` true while the table has no embedded
  `grhayl_nn_c2p` group. Route runtime replay coverage to Con2Prim tests and
  initialization behavior to EOS/Core initialization pages.

## Ownership Routes

- Core owns shared error-code vocabulary, `ghl_abort_if_error`, invalid routine
  name behavior, and `u0` helper routing. Start with
  [Core tests and fixtures](../core/tests-and-fixtures.md) and
  [Core errors, IO, debug, utilities](../core/errors-io-debug-utilities.md).
- EOS owns simple, hybrid, tabulated setup, HDF5 table reads, interpolation
  bounds, invalid EOS/table state, and beta-equilibrium root routes. Start with
  [EOS tests and fixtures](../gems/eos/tests-and-fixtures.md).
- Con2Prim owns invalid C2P key dispatch and recovery-selection semantics.
  Start with [Con2Prim tests and fixtures](../gems/con2prim/tests-and-fixtures.md)
  and [Con2Prim recovery flow](../gems/con2prim/recovery-flow.md).
- Con2Prim/EOS jointly route NN primitive-guess error keys: the loader lives in
  Con2Prim tabulated NN source, while key `85` is reached through tabulated EOS
  initialization.
- Neutrinos owns `Fermi` invalid-key coverage and NRPyLeakage fixture routes.
  Start with [Neutrinos tests and fixtures](../gems/neutrinos/tests-and-fixtures.md).

## Repo-Local References

- [configure](../../configure)
- [.github/run_tests.sh](../../.github/run_tests.sh)
- [.github/workflows/](../../.github/workflows/)
- [Unit_Tests/unit_test_code_error.c](../../Unit_Tests/unit_test_code_error.c)
- [Unit_Tests/nrpyleakage_main.h](../../Unit_Tests/nrpyleakage_main.h)
- [Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c](../../Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c)
- [Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c](../../Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c)
- [Unit_Tests/unit_test_nrpyleakage_luminosities.c](../../Unit_Tests/unit_test_nrpyleakage_luminosities.c)
- [Unit_Tests/unit_test_tabulated_eos.c](../../Unit_Tests/unit_test_tabulated_eos.c)
- [Unit_Tests/tabulated_eos_unit_test_helpers.c](../../Unit_Tests/tabulated_eos_unit_test_helpers.c)
- [Unit_Tests/test_compute_h_and_cs2.c](../../Unit_Tests/test_compute_h_and_cs2.c)
- [Unit_Tests/unit_test_tabulated_flux.c](../../Unit_Tests/unit_test_tabulated_flux.c)
- [Unit_Tests/unit_test_con2prim_tabulated.c](../../Unit_Tests/unit_test_con2prim_tabulated.c)
- [Unit_Tests/unit_test_con2prim_debug.c](../../Unit_Tests/unit_test_con2prim_debug.c)
- [Unit_Tests/unit_test_c2p_nn_guess.c](../../Unit_Tests/unit_test_c2p_nn_guess.c)
- [wiki/gems/eos/tests-and-fixtures.md](../gems/eos/tests-and-fixtures.md)
- [wiki/gems/con2prim/tests-and-fixtures.md](../gems/con2prim/tests-and-fixtures.md)
- [wiki/gems/neutrinos/tests-and-fixtures.md](../gems/neutrinos/tests-and-fixtures.md)
- [wiki/core/tests-and-fixtures.md](../core/tests-and-fixtures.md)
