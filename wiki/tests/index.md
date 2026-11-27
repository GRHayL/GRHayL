# Tests

> Router for `Unit_Tests` binaries, data generators, shared helpers, test runner outputs, reference data, HDF5 tables, ET_Legacy comparisons, error keys, and coverage gaps.

| Page | Go here when... |
| --- | --- |
| [Unit Test Catalog](unit-test-catalog.md) | You need the compact index of every `Unit_Tests/unit_test_*.c` and `Unit_Tests/compose/*.py` test with module links. |
| [Data Generators And Helpers](data-generators-and-helpers.md) | You need `Unit_Tests/data_gen/unit_test_data_*.c` generators, shared C/H helpers (`ghl_test_helpers.*`, `con2prim_newman_history.h`, `nrpyleakage_main.h`), or `pert_test_fail_*`. |
| [Runner And Build Outputs](runner-and-build-outputs.md) | You need `make tests`/`make datagen`, `.github/run_tests.sh`, `TestData` selector refs, `--disable-hdf5` filtering, or cleanup rules. |
| [Reference Data And Test Drivers](reference-data-and-test-drivers.md) | You need `ghl_pert_test_fail` comparison tolerances, reference-data source classes, and per-module binary data families. |
| [HDF5 Sample Tables](hdf5-sample-tables.md) | You need checked-in Hempel SFHo, generated `simple_table.h5`, downloaded LS220/SLy4 tables, `table_info.txt`, or `--disable-hdf5` table effects. |
| [ET Legacy Comparison](et-legacy-comparison.md) | You need `unit_test_ET_Legacy_*.c` replay semantics, limiter compatibility bounds, and the GRHayLib proof boundary. |
| [Expected Failure And Error Keys](expected-failure-and-error-keys.md) | You need `unit_test_code_error.c` keys `0`–`90`, inverted exit-code rules, and `GHL_DISABLE_HDF5` skip sets. |
| [Unit Test Coverage And Gap Matrix](unit-test-coverage-and-gap-matrix.md) | You need subsystem-by-subsystem coverage classification and documented test gaps across GRHayL. |

Back to [AGENTS.md](../../AGENTS.md).
