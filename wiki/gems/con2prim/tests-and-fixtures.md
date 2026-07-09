# Con2Prim Tests And Fixtures

This page routes Con2Prim-specific tests. For repo-wide test coverage, shared
test helpers, CI matrix notes, and non-Con2Prim fixtures, start with
`./wiki/test-map.md`.

Ground truth lives in `./Unit_Tests/`, `./Unit_Tests/data_gen/`,
`./.github/run_tests.sh`, and the workflow jobs that call the same test
binaries. Do not record fixture identity metadata here.

## Direct Tests

| Test | Behavior covered | Fixture/data dependency | Default runner status | HDF5/EOS-table needs | Relevant data generator |
| --- | --- | --- | --- | --- | --- |
| `./Unit_Tests/unit_test_con2prim_multi_method_hybrid.c` | Hybrid selected-method recovery for `Font1D`, `Palenzuela1D`, `Noble1D_entropy`, `Palenzuela1D_entropy`, `Noble1D`, and `Noble2D`; compares primitive outputs and return codes. | `metric_Bfield_initial_data.bin`, `con2prim_multi_method_hybrid_input.bin`, `con2prim_multi_method_hybrid_output.bin`, `con2prim_multi_method_hybrid_output_pert.bin`. `./.github/run_tests.sh` downloads these from `con2prim/` before running the binary. | Run by `./.github/run_tests.sh` as `./test/unit_test_con2prim_multi_method_hybrid`; also covered by workflow `c2p-routines` jobs. | No HDF5 table. Hybrid EOS setup is internal to the test. | `./Unit_Tests/data_gen/unit_test_data_con2prim_multi_method_hybrid.c` writes the shared metric/B-field and hybrid C2P fixtures. Ordinary runner path downloads fixtures; it does not regenerate them. |
| `./Unit_Tests/unit_test_con2prim_tabulated.c` | Tabulated multi-method recovery for `Palenzuela1D`, `Newman1D_entropy`, `Newman1D` with `Palenzuela1D` backup, `Palenzuela1D_entropy` with backup, and `Noble2D` with backup; checks main-vs-backup accounting and primitive tolerances across `rho_vs_T` and `Pmag_vs_Wm1` grids. With nonzero test key, it initializes `enable_neural_net_c2p` true, then replays the suite with neural-network primitive guesses disabled and enabled. | CLI EOS table path plus `con2prim_tabulated_*_{rho_vs_T,Pmag_vs_Wm1}_{unperturbed,perturbed}.bin` files for Palenzuela1D, Palenzuela1D_entropy, Newman1D, Newman1D_entropy, and Noble2D. `./.github/run_tests.sh` downloads these from `con2prim/` before running the binary. | Run by `./.github/run_tests.sh` as `./test/unit_test_con2prim_tabulated SLy4_3335_rho391_temp163_ye66.h5 1` after `pyghl append SLy4_3335_rho391_temp163_ye66.h5`; also covered by workflow `con-to-prim-tabulated` jobs with the same visible setup command. | Requires HDF5-enabled build and EOS table `SLy4_3335_rho391_temp163_ye66.h5`. The NN-enabled replay expects the table to contain an embedded `grhayl_nn_c2p` group after runner/workflow setup. `./configure --disable-hdf5` excludes tabulated tests. | Same test source has generation mode: `./test/unit_test_con2prim_tabulated <EOS table path> 0`. CI/default runner downloads fixtures, runs the visible `pyghl append` command, and invokes key `1`; do not document ordinary runs as regeneration or describe `pyghl` internals. |
| `./Unit_Tests/unit_test_c2p_nn_guess.c` | Direct tabulated NN primitive-guess coverage for model validation, bounded `x` inference fallback, root/embedded/legacy HDF5 model loads, malformed HDF5 dataset paths, and failed-load preservation of an existing model. | No external binary fixture. HDF5-enabled runs create temporary model files under `/tmp/unit_test_c2p_nn_*.h5`; no fixture download drives this direct test. | Run by `./.github/run_tests.sh` as `./test/unit_test_c2p_nn_guess`; also covered by workflow `c2p-failure` jobs through matrix entry `c2p_nn_guess`. | HDF5-enabled builds cover loader paths. With `GHL_DISABLE_HDF5`, the test still compiles the validation and inference-fallback checks while omitting HDF5 loader checks. | None. The HDF5 files are created by the test itself, not by a data generator. |
| `./Unit_Tests/unit_test_con2prim_debug.c` | Manual single-case tabulated C2P debug route using `ghl_con2prim_tabulated_multi_method` with `Newman1D`. Prints diagnostics/conservative/primitive state for investigation. | No binary fixture. Requires one CLI EOS table path. | Manual only. It is not called by `./.github/run_tests.sh`. It may be built in HDF5-enabled `make tests`, but no-HDF5 configure excludes it. | Requires HDF5-enabled build and caller-provided EOS table path. | None. |
| `./Unit_Tests/unit_test_hybrid_failure.c` | Expected hybrid C2P failure codes and backup-routine behavior, including `Noble2D` failures and `Font1D` backup success paths. | Hard-coded cases; no external binary fixture visible. | Run by `./.github/run_tests.sh` as `./test/unit_test_hybrid_failure`; also covered by workflow `c2p-failure` jobs. | No HDF5 table. Hybrid EOS setup is internal to the test. | None. |
| `./Unit_Tests/unit_test_apply_conservative_limits.c` | Conservative floor/limit application and diagnostic flags around `ghl_apply_conservative_limits`. | `metric_Bfield_initial_data.bin`, `apply_conservative_limits_input.bin`, `apply_conservative_limits_output.bin`, `apply_conservative_limits_output_pert.bin`. `./.github/run_tests.sh` downloads these from `con2prim/`. | Run by `./.github/run_tests.sh` as `./test/unit_test_apply_conservative_limits`; also covered by workflow `c2p-routines` jobs. | No HDF5 table. | `./Unit_Tests/data_gen/unit_test_data_con2prim_multi_method_hybrid.c` writes these fixtures. Ordinary runner path downloads fixtures; it does not regenerate them. |
| `./Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c` | Primitive floors/ceilings, entropy recomputation, speed limiting, and final `u0` calculation through `ghl_enforce_primitive_limits_and_compute_u0`. | `metric_Bfield_initial_data.bin`, `enforce_primitive_limits_and_compute_u0_input.bin`, `enforce_primitive_limits_and_compute_u0_output.bin`, `enforce_primitive_limits_and_compute_u0_output_pert.bin`. `./.github/run_tests.sh` downloads these from `con2prim/`. | Run by `./.github/run_tests.sh` as `./test/unit_test_enforce_primitive_limits_and_compute_u0`; also covered by workflow `c2p-routines` jobs. | No HDF5 table. | `./Unit_Tests/data_gen/unit_test_data_con2prim_multi_method_hybrid.c` writes these fixtures. Ordinary runner path downloads fixtures; it does not regenerate them. |
| `./Unit_Tests/unit_test_compute_conservs_and_Tmunu.c` | Primitive-to-conservative conversion and stress-energy tensor checks for `ghl_compute_conservs_and_Tmunu` and related conservative outputs. | `metric_Bfield_initial_data.bin`, `compute_conservs_and_Tmunu_input.bin`, `compute_conservs_and_Tmunu_output.bin`, `compute_conservs_and_Tmunu_output_pert.bin`. `./.github/run_tests.sh` downloads these from `con2prim/`. | Run by `./.github/run_tests.sh` as `./test/unit_test_compute_conservs_and_Tmunu`; also covered by workflow `c2p-routines` jobs. | No HDF5 table. | `./Unit_Tests/data_gen/unit_test_data_con2prim_multi_method_hybrid.c` writes these fixtures. Ordinary runner path downloads fixtures; it does not regenerate them. |
| `./Unit_Tests/unit_test_ET_Legacy_conservs.c` | Legacy comparison route for primitive-to-conservative data. Useful after conversion helper changes that may affect downstream compatibility. | `ET_Legacy_conservs_input.bin`, `ET_Legacy_conservs_output.bin`, `ET_Legacy_conservs_output_pert.bin`. `./.github/run_tests.sh` downloads these from `ET_Legacy/`. | Run by `./.github/run_tests.sh` as `./test/unit_test_ET_Legacy_conservs`; also covered by workflow `ET-Legacy` jobs. | No HDF5 table. | `./Unit_Tests/data_gen/unit_test_data_ET_Legacy_conservs.c` writes input-side files visible in source. Trusted output fixtures are downloaded by CI/default runner, not regenerated by ordinary test runs. |
| `./Unit_Tests/unit_test_ET_Legacy_primitives.c` | Legacy comparison route for conservative-to-primitive recovery. Useful after recovery-flow or solver behavior changes. | `ET_Legacy_primitives_input.bin`, `ET_Legacy_primitives_output.bin`, `ET_Legacy_primitives_output_pert.bin`. `./.github/run_tests.sh` downloads these from `ET_Legacy/`. | Run by `./.github/run_tests.sh` as `./test/unit_test_ET_Legacy_primitives`; also covered by workflow `ET-Legacy` jobs. | No HDF5 table. | `./Unit_Tests/data_gen/unit_test_data_ET_Legacy_primitives.c` writes input-side files visible in source. Trusted output fixtures are downloaded by CI/default runner, not regenerated by ordinary test runs. |

## Fixture Routes

- Shared hybrid/helper fixtures come from the `con2prim/` fixture namespace in
  `./.github/run_tests.sh`: `metric_Bfield_initial_data.bin`,
  `apply_conservative_limits_*`, `con2prim_multi_method_hybrid_*`,
  `enforce_primitive_limits_and_compute_u0_*`, and
  `compute_conservs_and_Tmunu_*`.
- Tabulated C2P fixtures also come from `con2prim/`, but they are paired with
  the EOS table downloaded as `SLy4_3335_rho391_temp163_ye66.h5.bz2` and
  decompressed before `unit_test_con2prim_tabulated` runs. The runner and
  workflows then run `pyghl append SLy4_3335_rho391_temp163_ye66.h5` before the
  NN-enabled tabulated replay; this page records only that visible setup
  command.
- `unit_test_con2prim_tabulated.c` has explicit generation mode key `0` and
  run mode key `1`. Current CI/default runner uses key `1` after downloading
  fixtures and after the visible `pyghl append` setup command.
- `unit_test_c2p_nn_guess.c` creates its own temporary HDF5 model files for
  direct loader tests. It is separate from the tabulated solver replay and does
  not consume downloaded `con2prim_tabulated_*` fixtures.
- `./configure --disable-hdf5` excludes `*tabulated*` unit tests,
  `unit_test_con2prim_debug.c`, and tabulated data generators from generated
  build lists.

## Targeted Checks After Build

Build first:

```sh
./configure -r
make tests
export LD_LIBRARY_PATH="$(pwd)/build/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
```

Use `./.github/run_tests.sh` as the full fixture-download-and-run route. For a
minimum local rerun, place the needed fixture files in repo root, matching the
filenames above, then run only the affected binaries:

| Change area | Minimum targeted route |
| --- | --- |
| Conservative limiting or `tau`/momentum floors | `./test/unit_test_apply_conservative_limits` |
| Hybrid solver dispatch, selected-method behavior, or hybrid method math | `./test/unit_test_con2prim_multi_method_hybrid` and `./test/unit_test_hybrid_failure` |
| Backup routing or failure diagnostics | `./test/unit_test_hybrid_failure`; add `./test/unit_test_con2prim_multi_method_hybrid` when success-path recovery may shift. |
| Primitive bounds, entropy recomputation, speed limit, or `u0` | `./test/unit_test_enforce_primitive_limits_and_compute_u0` |
| Primitive-to-conservative conversion or stress-energy outputs | `./test/unit_test_compute_conservs_and_Tmunu` and `./test/unit_test_ET_Legacy_conservs` |
| Conservative-to-primitive compatibility with legacy fixtures | `./test/unit_test_ET_Legacy_primitives` |
| Tabulated recovery, tabulated backup accounting, or HDF5 table routing | `./test/unit_test_con2prim_tabulated SLy4_3335_rho391_temp163_ye66.h5 1` |
| Direct NN primitive-guess validation, fallback, and HDF5 loader paths | `./test/unit_test_c2p_nn_guess` |
| One-off tabulated debug investigation | `./test/unit_test_con2prim_debug <EOS table path>`; manual only, not default runner coverage. |

When updating fixtures intentionally, use the matching data generator only after
reviewing its source and the expected fixture ownership. Do not infer that CI
regenerates fixtures from generator presence; current runner downloads them.
