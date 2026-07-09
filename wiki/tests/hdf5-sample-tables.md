# HDF5 Sample Tables

Purpose: route Unit_Tests HDF5 table assets, generated analytic tables,
downloaded runtime tables, no-HDF5 build effects, and downstream consumers.
Repo-local files remain authoritative; table API behavior belongs to the EOS
pages linked below.

Read with [EOS tests and fixtures](../gems/eos/tests-and-fixtures.md),
[EOS initialization and dispatch](../gems/eos/initialization-and-dispatch.md),
[tabulated table contract](../gems/eos/tabulated-table-contract.md),
[Con2Prim tests and fixtures](../gems/con2prim/tests-and-fixtures.md), and
[Neutrinos tests and fixtures](../gems/neutrinos/tests-and-fixtures.md).

## Table Roles

- Checked-in reduced Hempel table:
  [Unit_Tests/sample_table/Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817_simple.h5](../../Unit_Tests/sample_table/Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817_simple.h5)
  is a repo-local HDF5 sample table. It is reduced from a larger Hempel SFHo
  table and is useful as fixture evidence, not as the table used by the default
  full runner.
- Generated analytic table:
  [Unit_Tests/sample_table/generate_simple_table.py](../../Unit_Tests/sample_table/generate_simple_table.py)
  writes `simple_table.h5`. The script requires Python `numpy` and `h5py`,
  creates a small analytic table with dimensions `N_rho = 7`, `N_T = 5`, and
  `N_Y_e = 3`, and writes analytic datasets checked by
  `unit_test_tabulated_eos.c`.
- CI-downloaded analytic table:
  [.github/run_tests.sh](../../.github/run_tests.sh) downloads
  `EOS/simple_table.h5` from the GRHayL TestData route, places it in the repo
  root as `simple_table.h5`, and runs
  `./test/unit_test_tabulated_eos simple_table.h5`. Workflow files mirror this
  pattern with direct `curl` commands.
- Stellar-collapse LS220 table:
  `.github/run_tests.sh` downloads
  `LS220_234r_136t_50y_analmu_20091212_SVNr26.h5.bz2` from
  `https://stellarcollapse.org/EOS/`, unpacks
  `LS220_234r_136t_50y_analmu_20091212_SVNr26.h5`, then runs the tabulated
  Flux_Source test. `unit_test_tabulated_flux.c` hard-codes that local HDF5
  filename.
- Stellar-collapse SLy4 table:
  `.github/run_tests.sh` downloads
  `SLy4_3335_rho391_temp163_ye66.h5.bz2` from
  `https://stellarcollapse.org/EOS/`, unpacks
  `SLy4_3335_rho391_temp163_ye66.h5`, then uses that table for Neutrinos,
  Con2Prim tabulated recovery, and expected-error keys that need a table. The
  runner and workflows also run
  `pyghl append SLy4_3335_rho391_temp163_ye66.h5` before NN-enabled tabulated
  C2P replay; no repo-local page should describe `pyghl` internals or model
  provenance from that command alone.

## `table_info.txt`

[Unit_Tests/sample_table/table_info.txt](../../Unit_Tests/sample_table/table_info.txt)
records dimensions and ranges for the reduced Hempel sample, but it also
contains nonportable local absolute paths from the machine that created the
file. Do not copy those paths as authority. Use repo-local fixture filenames,
the generator script, runner downloads, workflow commands, and test sources
instead.

## Downstream Consumers

- `unit_test_tabulated_eos.c` consumes `simple_table.h5` and validates table
  dimensions, `energy_shift`, analytic table quantities, interpolation routes,
  `ghl_compute_h_and_cs2`, bounds, and beta-equilibrium helpers. Its analytic
  expectations are shared with `tabulated_eos_unit_test_helpers.c`.
- `test_compute_h_and_cs2.c` is helper-only code for tabulated tests. It
  assumes LS220-style table bounds and is not a standalone table contract.
- `unit_test_tabulated_flux.c` consumes
  `LS220_234r_136t_50y_analmu_20091212_SVNr26.h5` with downloaded
  `tabulated_flux_*` fixtures.
- `unit_test_con2prim_tabulated.c` consumes
  `SLy4_3335_rho391_temp163_ye66.h5` with downloaded
  `con2prim_tabulated_*` fixtures. Key `0` generates local data; the default
  runner uses key `1` after the visible `pyghl append` setup command, and the
  test replays both disabled and enabled NN primitive-guess modes when the key
  is nonzero.
- `unit_test_con2prim_debug.c` is a manual tabulated debug route that accepts a
  caller-provided EOS table path.
- `unit_test_nrpyleakage_*.c` use
  [Unit_Tests/nrpyleakage_main.h](../../Unit_Tests/nrpyleakage_main.h), which
  accepts an EOS table path plus key `0` or `1`. The default runner passes
  `SLy4_3335_rho391_temp163_ye66.h5` and key `1`.
- `unit_test_code_error.c` uses `SLy4_3335_rho391_temp163_ye66.h5` for
  HDF5-enabled tabulated setup, interpolation-bound, invalid table state, and
  root-bracketing error paths. NN key `85` uses the same table to cover the
  error path where NN initialization is enabled but the table lacks embedded
  `grhayl_nn_c2p`. Its table-read failure keys use `test.h5` instead: key `34`
  checks the missing-file path, while later read-table keys create malformed
  temporary `test.h5` inputs.
- `unit_test_c2p_nn_guess.c` creates temporary HDF5 model files under `/tmp`
  for root, embedded `grhayl_nn_c2p`, legacy, malformed-dataset, and
  failed-load-preservation cases. These are direct test artifacts, not external
  sample tables or downloaded binary fixtures.

## No-HDF5 Build Effects

No-HDF5 behavior is described only from [configure](../../configure). Passing
`./configure --disable-hdf5` adds `GHL_DISABLE_HDF5`, filters tabulated/HDF5
implementation sources, excludes `*tabulated*` unit tests, excludes
`unit_test_con2prim_debug.c`, excludes the three `unit_test_nrpyleakage_*.c`
tests, and filters tabulated data generators. `configure` says this disables
HDF5 and, for now, means no tabulated EOS.

Do not infer extra runtime exclusions from table filenames or workflow absence;
use `configure` for the build-list gate and use the test source for any
compiled expected-skip path.

## Repo-Local References

- [configure](../../configure)
- [.github/run_tests.sh](../../.github/run_tests.sh)
- [.github/workflows/](../../.github/workflows/)
- [Unit_Tests/sample_table/Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817_simple.h5](../../Unit_Tests/sample_table/Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817_simple.h5)
- [Unit_Tests/sample_table/generate_simple_table.py](../../Unit_Tests/sample_table/generate_simple_table.py)
- [Unit_Tests/sample_table/table_info.txt](../../Unit_Tests/sample_table/table_info.txt)
- [Unit_Tests/unit_test_tabulated_eos.c](../../Unit_Tests/unit_test_tabulated_eos.c)
- [Unit_Tests/tabulated_eos_unit_test_helpers.c](../../Unit_Tests/tabulated_eos_unit_test_helpers.c)
- [Unit_Tests/test_compute_h_and_cs2.c](../../Unit_Tests/test_compute_h_and_cs2.c)
- [Unit_Tests/unit_test_tabulated_flux.c](../../Unit_Tests/unit_test_tabulated_flux.c)
- [Unit_Tests/unit_test_con2prim_tabulated.c](../../Unit_Tests/unit_test_con2prim_tabulated.c)
- [Unit_Tests/unit_test_con2prim_debug.c](../../Unit_Tests/unit_test_con2prim_debug.c)
- [Unit_Tests/unit_test_c2p_nn_guess.c](../../Unit_Tests/unit_test_c2p_nn_guess.c)
- [Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c](../../Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c)
- [Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c](../../Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c)
- [Unit_Tests/unit_test_nrpyleakage_luminosities.c](../../Unit_Tests/unit_test_nrpyleakage_luminosities.c)
- [Unit_Tests/unit_test_code_error.c](../../Unit_Tests/unit_test_code_error.c)
- [Unit_Tests/nrpyleakage_main.h](../../Unit_Tests/nrpyleakage_main.h)
- [wiki/gems/eos/tests-and-fixtures.md](../gems/eos/tests-and-fixtures.md)
- [wiki/gems/con2prim/tests-and-fixtures.md](../gems/con2prim/tests-and-fixtures.md)
- [wiki/gems/neutrinos/tests-and-fixtures.md](../gems/neutrinos/tests-and-fixtures.md)
- [wiki/core/tests-and-fixtures.md](../core/tests-and-fixtures.md)
- Repo-local runner/workflow external table URLs:
  `https://raw.githubusercontent.com/GRHayL/TestData/main/EOS/simple_table.h5`,
  `https://stellarcollapse.org/EOS/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5.bz2`,
  and `https://stellarcollapse.org/EOS/SLy4_3335_rho391_temp163_ye66.h5.bz2`.
