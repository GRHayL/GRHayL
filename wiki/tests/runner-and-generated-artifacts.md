# Unit_Tests Runner And Generated Artifacts

This page routes how [Unit_Tests/](../../Unit_Tests) binaries are discovered, built, run,
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

[.github/run_tests.sh](../../.github/run_tests.sh) accepts
`[all|m1] [configure arguments...]` and passes the remaining arguments to
`./configure -r`. The default `all` suite builds with `make tests datagen`;
`m1` builds with `make tests`, compiling every configured test target but
executing only the Radiation M1 tests. Both modes export `LD_LIBRARY_PATH`
and `DYLD_LIBRARY_PATH` with `build/lib` before invoking compiled tests.
Most workflow jobs use the composite compile action, which runs
`make tests datagen` and `make install`; jobs then invoke binaries under
`test/`, often with `LD_LIBRARY_PATH` pointing at the installed `lib/`.
The `radiation-m1` jobs instead use the
[run_m1 action](../../.github/actions/run_m1/action.yml), which invokes the
runner's `m1` mode.

These stages are not interchangeable:

| Status | What repo evidence establishes |
| --- | --- |
| `targeted` | Generated Makefile selected source/binary. |
| `compiled` | Exact `make` child target completed and executable exists. |
| `compiled-unrun` | Binary exists, but no runner/workflow command executed it. |
| `executed` | Exact binary ran with recorded arguments and prerequisites. |
| `workflow-only` | YAML contains a job/matrix command; presence gives no local execution result. |

Aggregate `make tests`/`make datagen` success is compile/link evidence only.
`make datagen` never generates fixture data unless its executables are then run.

## Artifact Classes

| Class | Normal location | Source evidence | Notes |
| --- | --- | --- | --- |
| Build products | `build/`, `build/lib/`, `test/` | `configure` generated Makefile rules | `test/` contains compiled unit-test executables from `make tests`. |
| Data-generator executables | `test/data_gen/` | `configure` generated Makefile rules | Built by `make datagen`; not run by default by [.github/run_tests.sh](../../.github/run_tests.sh). |
| Downloaded root-level fixtures | repo root, usually `*.bin`, `*.h5`, `*.bz2` before decompression | `.github/run_tests.sh` and workflows | Runner downloads binary fixtures from the repo-visible `GRHayL/TestData` paths and EOS tables from visible `stellarcollapse.org/EOS` URLs. |
| Generated root-level fixtures | repo root, names chosen by generator/test source | `Unit_Tests/data_gen/*`, [Unit_Tests/nrpyleakage_main.h](../../Unit_Tests/nrpyleakage_main.h), [Unit_Tests/unit_test_con2prim_tabulated.c](../../Unit_Tests/unit_test_con2prim_tabulated.c) | Data generators and generation-mode tests write local fixture files when invoked manually. |
| Checked-in sample inputs | [Unit_Tests/sample_table/](../../Unit_Tests/sample_table) | [wiki/test-map.md](../test-map.md), [wiki/generated-boundaries.md](../generated-boundaries.md), test sources | The reduced HDF5 sample table and table generator are repo files, not outputs from `configure`. |

Examples of generation-mode tests from repo-local test sources:
[Unit_Tests/nrpyleakage_main.h](../../Unit_Tests/nrpyleakage_main.h) uses key
`0` for data generation and key `1` for unit-test replay, while
[Unit_Tests/unit_test_con2prim_tabulated.c](../../Unit_Tests/unit_test_con2prim_tabulated.c)
generates `con2prim_tabulated_*_{unperturbed,perturbed}.bin` files when run in
its generation mode. Normal CI replay downloads fixtures and runs test mode.

## Runtime Runner

[.github/run_tests.sh](../../.github/run_tests.sh) is the broad local replay
route, not a complete workflow matrix. The default `all` suite follows this
sequence:

1. Run `./configure -r`.
2. Run `make tests datagen`.
3. Export `LD_LIBRARY_PATH` and `DYLD_LIBRARY_PATH` with `build/lib`.
4. Run the `unit_test_coverage_*` boundary executables selected by the runner.
5. Download missing root-level `*.bin` fixtures and HDF5/EOS tables, reusing
   existing regular files.
6. Decompress downloaded `*.bz2` EOS tables.
7. Run selected binaries under `test/`, including the direct
   `unit_test_c2p_nn_guess` route.
8. Run `unit_test_code_error` keys `0` through `90` as expected process
   failures.
9. Run `pyghl append SLy4_3335_rho391_temp163_ye66.h5` before the tabulated
   C2P replay. This page records only that visible setup command, not `pyghl`
   internals.
10. Continue selected compiled-test runs, then call `run_m1_tests`.
11. Use an `EXIT` trap to remove tracked download/decompression paths, the
    expected-error work directory, and runner-owned M1 temporary directories.

The `m1` suite runs `make tests` and only the Radiation M1 executables, then
exits before legacy fixture downloads, expected-error-key runs, and `pyghl`
setup. With `M1_FIXTURE_DIR` unset, the runner downloads every named
`radiation/*.bin.gz` member at the pinned revision recorded in
[.github/radiation-testdata-ref](../../.github/radiation-testdata-ref) and
prepares them in a private temporary directory. Supplying `M1_FIXTURE_DIR`
with a raw or gzip
package bypasses downloads; a complete raw `.bin` directory is used directly,
otherwise [prepare_m1_fixtures.sh](../../.github/prepare_m1_fixtures.sh)
prepares members in a runner-owned temporary directory. Direct executable
invocations without fixture arguments remain local-only. The provider test
uses `--generated-fixture` for its local input. See
[M1 tests and fixtures](../gems/radiation-m1/tests-and-fixtures.md) for replay
arguments and fixture ownership. These fixture members are published
`GRHayL/TestData` payloads described by that repository's provenance route;
replay consumes them and never regenerates expected values locally.

The [Radiation CI action](../../.github/actions/run_m1/action.yml) forwards
its optional `fixture-dir` input as `M1_FIXTURE_DIR` and invokes the runner's
`m1` mode. It performs no TestData checkout; the runner retries the pinned
download by default. The pinned revision must be available from TestData.

The cleanup trap runs on success or failure and deletes only its tracked
paths, including runner-owned M1 download and prepared-fixture temporary
directories. Existing regular fixture files (or EOS archives/uncompressed files)
bypass downloads, so the configured URL does not establish their provenance.
This cleanup ownership does not promise unchanged contents: configure/build
steps write artifacts, and `pyghl append` receives the SLy4 table even when it
preexisted. The optically thin gas test also creates, truncates, and removes
its fixed marker filename in the current directory. Error key `85` requires that table to lack `grhayl_nn_c2p` before
the later append step. Use a disposable checkout for the full runner, as in
[Workflows](../workflows.md); arbitrary path collisions are not a supported
preservation mechanism.
`make realclean` separately removes build products such as `build/`, `test/`,
and the generated `Makefile`.

`make clean` is narrower: generated rules omit individual unit-test objects
(`TEXEOBJS`), data-generator executables, and data-generator objects from its
removal list. `make realclean` removes
their parent `test/` and build directories. In HDF5 mode,
`unit_test_c2p_nn_guess` owns a unique private temporary directory and enumerates and unlinks its contents
before removing the directory on normal or handled failure exit.

## Fixture Reference Selection

The runner reads [.github/con2prim-testdata-ref](../../.github/con2prim-testdata-ref)
for shared hybrid Con2Prim/Core fixtures and tabulated Con2Prim replay fixtures.
Its other TestData downloads use
[.github/et-legacy-testdata-ref](../../.github/et-legacy-testdata-ref), including
ET Legacy, analytic EOS, Core, Flux_Source, PLM, Neutrinos, and Induction data.
The M1 replay route reads the separate
[.github/radiation-testdata-ref](../../.github/radiation-testdata-ref)
selector for the published `radiation/*.bin.gz` members.
Link these selector files rather than copying their values into the KB.

Workflow download selection differs: shared hybrid Con2Prim/Core fixtures use
the Con2Prim selector, and ET Legacy/Flux_Source use the ET Legacy selector.
Standalone Core, analytic EOS, Reconstruction, Induction, Neutrinos, and
tabulated Con2Prim jobs download from TestData `main`. Inspect the owning
[workflow](../../.github/workflows/) for each command. Neither a download
command nor a current local generator proves how a historical fixture was
produced; preexisting-file reuse adds another local provenance boundary.

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

Treat `.github/run_tests.sh` as one broad local-style driver, not a complete
enumeration of every workflow job. Treat workflows as CI matrices, not proof
that normal local runs regenerate trusted fixtures.

Set comparison is exact in current sources: default `configure` targets the
`Unit_Tests/unit_test_*.c` set, while the runner directly invokes all except
`unit_test_WENOZ_reconstruction`,
`unit_test_tabulated_eos_compose`, and `unit_test_con2prim_debug` — the `all`
suite finishes with `run_m1_tests`, so the Radiation M1 executables are part
of the aggregate default replay; workflows
invoke WENOZ and the focused CompOSE test, while no normal runner/workflow
invocation for the debug binary is visible.

NN primitive-guess coverage appears in both paths: `.github/run_tests.sh` runs
`./test/unit_test_c2p_nn_guess`, while workflow `c2p-failure` matrices include
`c2p_nn_guess`. The tabulated replay path is separate: runner and workflows run
`pyghl append SLy4_3335_rho391_temp163_ye66.h5` before
`./test/unit_test_con2prim_tabulated SLy4_3335_rho391_temp163_ye66.h5 1`.

## HDF5 Filtering

By default, `configure` expects HDF5 and includes tabulated/HDF5 tests and data
generators. With `--disable-hdf5`, `configure` adds `-DGHL_DISABLE_HDF5`,
applies its exact implementation-source predicate, and changes generated test
and data-generator target lists. That predicate retains
`tabulated_primitive_guess_helpers.c` and the
`Con2Prim/Tabulated/neural_network_guess/` sources, plus the concrete tabulated
and tabulated-entropy HLLE flux kernels despite their paths. The test filter:

- Excludes `unit_test_*tabulated*.c`.
- Excludes `unit_test_con2prim_debug.c`.
- Excludes `unit_test_nrpyleakage_constant_density_sphere.c`.
- Excludes `unit_test_nrpyleakage_luminosities.c`.
- Excludes `unit_test_nrpyleakage_optically_thin_gas.c`.
- Excludes `Unit_Tests/data_gen/*tabulated*`.

The table-free NRPyLeakage physics and classifier-fallback tests remain build
targets, as does `unit_test_c2p_nn_guess` with its HDF5 portions compiled out.
The focused no-HDF5 `c2p_nn_guess` leg in
[Ubuntu-Clang](../../.github/workflows/github-actions-Ubuntu-clang.yml) runs
`make tests install`, inspects the direct tabulated HLLE exported definitions,
and executes the NN test. This does not execute every compiled test or prove
an external consumer can link those definitions. The aggregate
`.github/run_tests.sh` uses default `./configure -r`; it is not a no-HDF5 driver.

[Unit_Tests/unit_test_code_error.c](../../Unit_Tests/unit_test_code_error.c) has a separate no-HDF5 behavior: when
compiled with `GHL_DISABLE_HDF5`, HDF5-only error keys are treated as skipped
passes. That is test-source behavior; the target filtering above is from
`configure`.

## Test-Specific Evidence For This Route

These test files determine the classification of generation modes,
manual-only tests, HDF5 needs, and expected failures:

- [Unit_Tests/nrpyleakage_main.h](../../Unit_Tests/nrpyleakage_main.h)
- [Unit_Tests/unit_test_con2prim_tabulated.c](../../Unit_Tests/unit_test_con2prim_tabulated.c)
- [Unit_Tests/unit_test_con2prim_debug.c](../../Unit_Tests/unit_test_con2prim_debug.c)
- [Unit_Tests/unit_test_c2p_nn_guess.c](../../Unit_Tests/unit_test_c2p_nn_guess.c)
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
