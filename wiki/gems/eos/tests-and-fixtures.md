# EOS Tests And Fixtures

## Purpose

This page routes EOS-specific tests, helper-only files, sample tables, generated
fixtures, CI downloads, and error-path coverage. Source, tests, CI scripts, and
build configuration remain authoritative; this page separates direct EOS
coverage from dependent impact signals.

Read with [EOS initialization and dispatch](initialization-and-dispatch.md),
[hybrid piecewise-polytrope](hybrid-piecewise-polytrope.md),
[tabulated table contract](tabulated-table-contract.md), and
[tabulated interpolation and bounds](tabulated-interpolation-and-bounds.md).

## Direct EOS Tests

- [Unit_Tests/unit_test_piecewise_polytrope.c](../../../Unit_Tests/unit_test_piecewise_polytrope.c)
  directly initializes a hybrid piecewise-polytrope EOS and checks derived
  `K_ppoly` and `eps_integ_const` entries. This is the direct regression route
  for piecewise-polytrope constant setup. It does not check the return code,
  `p_ppoly`, or invalid-input behavior; see the hybrid contract's breakpoint
  contradiction.
- [Unit_Tests/unit_test_tabulated_eos.c](../../../Unit_Tests/unit_test_tabulated_eos.c)
  directly initializes a tabulated EOS from a CLI HDF5 table path, checks table
  dimensions and `energy_shift`, validates analytic table quantities,
  exercises direct `(rho,Y_e,T)` wrappers, auxiliary wrappers from `eps`, `P`,
  and `S`, `ghl_compute_h_and_cs2`, selected bound-enforcement paths, and
  beta-equilibrium rho-map helpers.
- [Unit_Tests/unit_test_tabulated_eos_compose.c](../../../Unit_Tests/unit_test_tabulated_eos_compose.c)
  independently reads a converter-produced StellarCollapse file, checks every
  node of all 19 serialized fields against GRHayL's loader mapping and units,
  validates derived enthalpy and the relativistic sound-speed route, exercises
  midpoint storage-space interpolation and `eps/P/S/h` inversions, checks the
  six-value NRPyLeakage interpolation order and range errors, requires a
  no-fallback Palenzuela recovery, compares characteristic speeds, tabulated
  HLLE and entropy fluxes, and source terms with analytic goldens, compares all
  eight combined-leakage outputs with an independent high-precision
  Ruffert-equation oracle, and frees table memory. `configure` auto-discovers
  it only in HDF5-enabled builds.
- [Unit_Tests/unit_test_code_error.c](../../../Unit_Tests/unit_test_code_error.c)
  covers EOS initialization errors for simple, hybrid, and tabulated setup;
  tabulated interpolation/helper errors including out-of-table and
  too-many-variable cases; table-read errors; invalid table/EOS state; and
  beta-equilibrium rho-map error paths.

The tabulated test directly calls many, but not all, registered wrappers.
[Tabulated interpolator catalog](tabulated-interpolator-catalog.md) separates
declaration, storage, assignment, build, and direct-call evidence.

The Python suite under
[Unit_Tests/compose/](../../../Unit_Tests/compose/) uses an asymmetric
official-schema-shaped CompOSE fixture. It covers fixed profile and control
validation, HDF5 structure and security boundaries, every output mapping and
regularization branch, publication failures and cleanup, and every CLI result.
The focused Ubuntu GCC workflow requires 100% production Python line and
branch coverage and uses that fixture to feed the C integration test. The full
table-141 artifact remains external and is a manual qualification input, not a
repository fixture.

## Helper-Only Files

- [Unit_Tests/test_compute_h_and_cs2.c](../../../Unit_Tests/test_compute_h_and_cs2.c)
  provides a test-only `ghl_test_compute_h_and_cs2` helper for enthalpy and
  sound-speed checks. It assumes LS220-style table bounds and is not a
  standalone unit test.
- [Unit_Tests/tabulated_eos_unit_test_helpers.c](../../../Unit_Tests/tabulated_eos_unit_test_helpers.c)
  provides analytic table-quantity helpers used by tabulated EOS validation. It
  is helper code, not a standalone test binary.

## Error And HDF5 Behavior

`unit_test_code_error.c` has two HDF5 modes:

- HDF5-enabled builds run HDF5/tabulated error keys, including tabulated EOS
  setup failures, direct and auxiliary interpolation bound errors,
  too-many-variable helper calls, table-read failures, invalid EOS/table type,
  and beta-equilibrium rho-map bounds.
- No-HDF5 builds compile with `GHL_DISABLE_HDF5`; the test treats HDF5-only
  keys as skipped passes. Configuration excludes most tabulated/HDF5 sources
  and tabulated-named tests, while retaining Con2Prim NN validation/inference
  sources and `unit_test_c2p_nn_guess`; its HDF5 loader cases are compiled out.

The split is source-backed by
[Unit_Tests/unit_test_code_error.c](../../../Unit_Tests/unit_test_code_error.c),
[configure](../../../configure), [README.md](../../../README.md), and
[wiki/build-and-ci.md](../../build-and-ci.md).

## Table And Fixture Roles

- Checked-in reduced sample table:
  [Unit_Tests/sample_table/Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817_simple.h5](../../../Unit_Tests/sample_table/Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817_simple.h5)
  is a repo-local HDF5 sample table. Its dimensions/ranges are described by
  [Unit_Tests/sample_table/table_info.txt](../../../Unit_Tests/sample_table/table_info.txt).
- Generated analytic table:
  [Unit_Tests/sample_table/generate_simple_table.py](../../../Unit_Tests/sample_table/generate_simple_table.py)
  writes `simple_table.h5` with analytic quantities used by
  `unit_test_tabulated_eos.c` comparisons. The generated file name differs from
  the checked-in Hempel sample file.
- CI/runtime downloaded table:
  [.github/run_tests.sh](../../../.github/run_tests.sh) downloads
  `EOS/simple_table.h5` from the GRHayL test-data route before running
  `./test/unit_test_tabulated_eos simple_table.h5`. The same runner downloads
  larger stellar-collapse tables for dependent tabulated Flux_Source,
  Con2Prim, and Neutrinos tests.

`table_info.txt` records local absolute paths from the environment that created
the reduced Hempel sample. Do not copy those paths as authority; use repo-local
fixture files, generator source, CI runner paths, and table-reading code
instead.

## Dependent Impact Tests

Dependent tests are impact signals for EOS changes, not primary EOS contracts:

- [Unit_Tests/unit_test_grhayl_core_test_suite.c](../../../Unit_Tests/unit_test_grhayl_core_test_suite.c)
  covers `ghl_set_prims_to_constant_atm` for simple and hybrid EOS setup. Its
  tabulated Atmosphere branch is commented out, and the file has a TODO to add
  table-based default checks before extending the loop. For simple/hybrid it
  does not assert `Y_e` or temperature, even though constant Atmosphere copies
  `eos.Y_e_atm` and `eos.T_atm` and those initializers do not set the fields.
  Its simple-EOS default-floor check also does not assert `eps_min` or
  `entropy_min` after both density and pressure minima default to zero.
- Con2Prim tabulated tests route through
  [Con2Prim tests and fixtures](../con2prim/tests-and-fixtures.md).
- Neutrinos table-backed tests route through
  [Neutrinos tests and fixtures](../neutrinos/tests-and-fixtures.md).
- Tabulated Flux_Source tests route through [Flux Source](../flux-source.md)
  and [Test Map](../../test-map.md).

## Source Of Truth

- [Unit_Tests/unit_test_piecewise_polytrope.c](../../../Unit_Tests/unit_test_piecewise_polytrope.c)
- [Unit_Tests/unit_test_tabulated_eos.c](../../../Unit_Tests/unit_test_tabulated_eos.c)
- [Unit_Tests/unit_test_tabulated_eos_compose.c](../../../Unit_Tests/unit_test_tabulated_eos_compose.c)
- [Unit_Tests/compose/](../../../Unit_Tests/compose/)
- [tools/compose/compose_to_grhayl.py](../../../tools/compose/compose_to_grhayl.py)
- [Unit_Tests/test_compute_h_and_cs2.c](../../../Unit_Tests/test_compute_h_and_cs2.c)
- [Unit_Tests/tabulated_eos_unit_test_helpers.c](../../../Unit_Tests/tabulated_eos_unit_test_helpers.c)
- [Unit_Tests/unit_test_code_error.c](../../../Unit_Tests/unit_test_code_error.c)
- [Unit_Tests/unit_test_grhayl_core_test_suite.c](../../../Unit_Tests/unit_test_grhayl_core_test_suite.c)
- [Unit_Tests/sample_table/Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817_simple.h5](../../../Unit_Tests/sample_table/Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817_simple.h5)
- [Unit_Tests/sample_table/generate_simple_table.py](../../../Unit_Tests/sample_table/generate_simple_table.py)
- [Unit_Tests/sample_table/table_info.txt](../../../Unit_Tests/sample_table/table_info.txt)
- [.github/run_tests.sh](../../../.github/run_tests.sh)
- [.github/workflows/](../../../.github/workflows/)
- [README.md](../../../README.md)
- [configure](../../../configure)
- [wiki/test-map.md](../../test-map.md)
- [wiki/build-and-ci.md](../../build-and-ci.md)
