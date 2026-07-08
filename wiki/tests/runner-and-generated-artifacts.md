# Unit_Tests Runner And Generated Artifacts

This page routes how `Unit_Tests/` binaries are discovered, built, run,
HDF5-filtered, and cleaned. Repo-local runner/build evidence is authoritative:
[configure](../../configure), [.github/run_tests.sh](../../.github/run_tests.sh),
and [.github/workflows](../../.github/workflows/). Unit-test sources explain
test-specific modes and file names; runner behavior comes from the scripts and
workflow files.

## Configure And Build Targets

[configure](../../configure) discovers unit-test targets from
`Unit_Tests/unit_test_*.c` and data-generator targets from
`Unit_Tests/data_gen/unit_test_data_*.c`. It writes a generated `Makefile` with
`tests` and `datagen` targets:

```sh
./configure -r
make tests
make datagen
```

`make tests` builds test executables under `test/`. `make datagen` builds data
generator executables under `test/data_gen/`. Both link against the configured
GRHayL shared library under `build/lib/`, so local direct runs need a runtime
library path like:

```sh
export LD_LIBRARY_PATH="$(pwd)/build/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
```

[.github/run_tests.sh](../../.github/run_tests.sh) runs `./configure -r`, then
`make tests`, then exports `LD_LIBRARY_PATH` with `build/lib` before invoking
compiled tests. Workflow jobs use the composite compile action, which runs
`make tests datagen` and `make install`; jobs then invoke binaries under
`test/`, often with `LD_LIBRARY_PATH` pointing at the installed `lib/`.

## Artifact Classes

| Class | Normal location | Source evidence | Notes |
| --- | --- | --- | --- |
| Build products | `build/`, `build/lib/`, `test/` | `configure` generated Makefile rules | `test/` contains compiled unit-test executables from `make tests`. |
| Data-generator executables | `test/data_gen/` | `configure` generated Makefile rules | Built by `make datagen`; not run by default by `.github/run_tests.sh`. |
| Downloaded root-level fixtures | repo root, usually `*.bin`, `*.h5`, `*.bz2` before decompression | `.github/run_tests.sh` and workflows | Runner downloads binary fixtures from the repo-visible `GRHayL/TestData` paths and EOS tables from visible `stellarcollapse.org/EOS` URLs. |
| Generated root-level fixtures | repo root, names chosen by generator/test source | `Unit_Tests/data_gen/*`, `Unit_Tests/nrpyleakage_main.h`, `Unit_Tests/unit_test_con2prim_tabulated.c` | Data generators and generation-mode tests write local fixture files when invoked manually. |
| Checked-in sample inputs | `Unit_Tests/sample_table/` | `wiki/test-map.md`, `wiki/generated-boundaries.md`, test sources | The reduced HDF5 sample table and table generator are repo files, not outputs from `configure`. |

Examples of generation-mode tests from repo-local test sources:
[Unit_Tests/nrpyleakage_main.h](../../Unit_Tests/nrpyleakage_main.h) uses key
`0` for data generation and key `1` for unit-test replay, while
[Unit_Tests/unit_test_con2prim_tabulated.c](../../Unit_Tests/unit_test_con2prim_tabulated.c)
generates `con2prim_tabulated_*_{unperturbed,perturbed}.bin` files when run in
its generation mode. Normal CI replay downloads fixtures and runs test mode.

## Runtime Runner

[.github/run_tests.sh](../../.github/run_tests.sh) is the compact full-suite
local route:

1. Run `./configure -r`.
2. Run `make tests`.
3. Export `LD_LIBRARY_PATH` with `build/lib`.
4. Download root-level `*.bin` fixtures and HDF5/EOS tables.
5. Decompress downloaded `*.bz2` EOS tables.
6. Run selected binaries under `test/`.
7. Run `unit_test_code_error` keys `0` through `82` as expected process
   failures.
8. Remove root-level `*.bin`, `*.h5`, and `*.bz2`.

The cleanup is explicit:

```sh
rm -f ./*.bin ./*.h5 ./*.bz2
```

If the runner exits early, root-level downloaded or decompressed fixtures may
remain. `make realclean` removes build products such as `build/`, `test/`, and
the generated `Makefile`; it does not own remote fixture cleanup.

## Workflow Matrix Differences

Workflow matrices can run tests absent from `.github/run_tests.sh`. The
clearest visible case is WENOZ reconstruction: `.github/run_tests.sh` downloads
and runs `PLM_reconstruction` but does not call `unit_test_WENOZ_reconstruction`;
the workflow `reconstruction` matrices include both `PLM_reconstruction` and
`WENOZ_reconstruction`. Source confirmation lives in
[Unit_Tests/unit_test_WENOZ_reconstruction.c](../../Unit_Tests/unit_test_WENOZ_reconstruction.c)
and [Unit_Tests/data_gen/unit_test_data_WENOZ_reconstruction.c](../../Unit_Tests/data_gen/unit_test_data_WENOZ_reconstruction.c),
with routing detail in
[Reconstruction tests and fixtures](../gems/reconstruction/tests-and-fixtures.md).

Treat `.github/run_tests.sh` as one full local-style driver, not a complete
enumeration of every workflow job. Treat workflows as CI matrices, not proof
that normal local runs regenerate trusted fixtures.

## HDF5 Filtering

By default, `configure` expects HDF5 and includes tabulated/HDF5 tests and data
generators. With `--disable-hdf5`, `configure` adds `-DGHL_DISABLE_HDF5`,
filters tabulated implementation sources, and changes the generated test and
data-generator target lists:

- Excludes `unit_test_*tabulated*.c`.
- Excludes `unit_test_con2prim_debug.c`.
- Excludes `unit_test_nrpyleakage_constant_density_sphere.c`.
- Excludes `unit_test_nrpyleakage_luminosities.c`.
- Excludes `unit_test_nrpyleakage_optically_thin_gas.c`.
- Excludes `Unit_Tests/data_gen/*tabulated*`.

`Unit_Tests/unit_test_code_error.c` has a separate no-HDF5 behavior: when
compiled with `GHL_DISABLE_HDF5`, HDF5-only error keys are treated as skipped
passes. That is test-source behavior; the target filtering above is from
`configure`.

## Test-Specific Evidence Read For This Route

These test files were checked to avoid misclassifying generation modes,
manual-only tests, HDF5 needs, and expected failures:

- [Unit_Tests/nrpyleakage_main.h](../../Unit_Tests/nrpyleakage_main.h)
- [Unit_Tests/unit_test_con2prim_tabulated.c](../../Unit_Tests/unit_test_con2prim_tabulated.c)
- [Unit_Tests/unit_test_con2prim_debug.c](../../Unit_Tests/unit_test_con2prim_debug.c)
- [Unit_Tests/unit_test_tabulated_eos.c](../../Unit_Tests/unit_test_tabulated_eos.c)
- [Unit_Tests/unit_test_code_error.c](../../Unit_Tests/unit_test_code_error.c)

## Repo-Local References

- [AGENTS.md](../../AGENTS.md)
- [wiki/index.md](../index.md)
- [wiki/catalog.md](../catalog.md)
- [wiki/test-map.md](../test-map.md)
- [wiki/build-and-ci.md](../build-and-ci.md)
- [wiki/generated-boundaries.md](../generated-boundaries.md)
- [configure](../../configure)
- [.github/run_tests.sh](../../.github/run_tests.sh)
- [.github/workflows](../../.github/workflows/)
- [.github/actions/compile_GRHayL/action.yml](../../.github/actions/compile_GRHayL/action.yml)

## External References

None used. Repo-local files above are the authority for this page.
