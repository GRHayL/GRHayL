# Runner And Build Outputs

> How `configure`, `make tests`/`make datagen`, `.github/run_tests.sh`, and CI workflows select, download, run, filter, and clean `Unit_Tests` binaries and data. · Status: confirmed
> Up: [Tests](index.md)

## Summary

`configure` generates `tests` (`test/`) and `datagen` (`test/data_gen/`) targets from `Unit_Tests/unit_test_*.c` and `Unit_Tests/data_gen/unit_test_data_*.c`, filtered by `--disable-hdf5`. `.github/run_tests.sh` runs a broad local replay pass (downloading missing `GRHayL/TestData` files and StellarCollapse tables into repo root and cleaning tracked downloads on `EXIT`), but omits WENO-z reconstruction, the CompOSE integration test, and manual Con2Prim debug binaries.

## Detail

### Build Targets, Output Classes, And Execution Stages

`make tests` and `make datagen` compile and link binaries against `<builddir>/lib/libghl*` (`<builddir>` defaults to `build`; see [Build And Configure](../architecture/build-and-configure.md#configure-options-and-floating-point-guards)); neither executes tests or generates reference data. To run them locally on Linux, prepend the absolute path of `<builddir>/lib` to `LD_LIBRARY_PATH` (for example `LD_LIBRARY_PATH="$(pwd)/build/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"` after `./configure -r`). Distinguish `targeted` (in generated `Makefile`), `compiled`, `compiled-unrun`, `executed`, and `workflow-only` evidence.

| Output class | Location | Cleanup & provenance boundary |
| --- | --- | --- |
| Library & test binaries | `<builddir>/`, `<builddir>/lib/`, `test/`, `test/data_gen/` | `make clean` keeps `TEXEOBJS` and `data_gen` binaries/objects; `make realclean` deletes `<builddir>/` (default `build/`), `test/`, and `Makefile`. (Ev: `configure::realclean`) |
| Downloaded root data | Repo root (`*.bin`, `*.h5`, `*.bz2`) | `.github/run_tests.sh` tracks only files it downloads/decompresses plus its private `code_error_workdir` for `EXIT` cleanup. (Ev: `run_tests.sh::code_error_workdir`) |
| Generation-mode outputs | Repo root (`*_unperturbed.bin`, `*_perturbed.bin`) | Written only when generators or key-`0` modes (`nrpyleakage_main.h`, `unit_test_con2prim_tabulated.c`) are invoked manually. (Ev: `nrpyleakage_main.h`) |
| Checked-in sample data | `Unit_Tests/sample_table/` | Tracked repository inputs, not `configure` outputs. (Ev: `generate_simple_table.py`) |

### `.github/run_tests.sh`, Reference Selectors, And Matrix Differences

Claim status: contested; contradiction: CONTR-0003. See [CONTR-0003](../contradictions.md#contr-0003).

`.github/run_tests.sh` executes `./configure -r`, `make tests datagen`, exports `LD_LIBRARY_PATH`, downloads missing root-level `*.bin` and EOS `.h5.bz2` files (reusing any preexisting regular files without checking provenance), decompresses `.bz2` archives, runs compiled tests including `unit_test_c2p_nn_guess`, runs `unit_test_code_error` over keys `0`–`90` inside its private `code_error_workdir`, runs `pyghl append SLy4_3335_rho391_temp163_ye66.h5` before `./test/unit_test_con2prim_tabulated SLy4_3335_rho391_temp163_ye66.h5 1`, and cleans tracked downloads on `EXIT`. Always run `.github/run_tests.sh` inside a disposable checkout: `pyghl append` mutates a reused SLy4 table (while key `85` requires an unmodified SLy4 table without `grhayl_nn_c2p`), and `unit_test_nrpyleakage_optically_thin_gas` writes and removes a fixed marker filename in cwd. In HDF5 mode, `unit_test_c2p_nn_guess` manages and removes its own private temporary directory.

- **TestData reference selectors:** `.github/run_tests.sh` reads `.github/con2prim-testdata-ref` for shared hybrid Con2Prim/Core and tabulated Con2Prim data, and `.github/et-legacy-testdata-ref` for ET_Legacy, analytic EOS, Core, Flux_Source, PLM, Neutrinos, and Induction data. Workflows use those selectors for Con2Prim/ET_Legacy/Flux_Source while downloading standalone Core, analytic EOS, Reconstruction, Induction, and Neutrinos files from `TestData` `main`.
- **Omissions from `run_tests.sh`:** Default `make tests` builds every `Unit_Tests/unit_test_*.c`, whereas `.github/run_tests.sh` skips `unit_test_WENOZ_reconstruction` (run in CI `reconstruction` jobs), `unit_test_tabulated_eos_compose` (run in Ubuntu-GCC `compose-regularized-eos`), and `unit_test_con2prim_debug` (manual-only).

Evidence (CI behavior): `run_tests.sh::con2prim-testdata-ref`; corroboration: `et-legacy-testdata-ref`.

### `--disable-hdf5` Test Filtering

Passing `./configure --disable-hdf5` adds `-DGHL_DISABLE_HDF5`, retains `tabulated_primitive_guess_helpers.c`, `Con2Prim/Tabulated/neural_network_guess/`, disabled direct-C2P stubs, and concrete tabulated HLLE flux kernels, and excludes `unit_test_*tabulated*.c`, `unit_test_con2prim_debug.c`, `unit_test_nrpyleakage_{constant_density_sphere,luminosities,optically_thin_gas}.c`, and `Unit_Tests/data_gen/*tabulated*`. Table-free NRPyLeakage `physics`/`classifier_fallback` and `unit_test_c2p_nn_guess` (with HDF5 paths compiled out) remain built; `unit_test_code_error.c` skips HDF5-only keys internally. `.github/run_tests.sh` uses default HDF5 `./configure -r` and is not a no-HDF5 driver.

Evidence (descriptive behavior): `configure::GHL_DISABLE_HDF5`; corroboration: `github-actions-Ubuntu-clang.yml::c2p_nn_guess`.

## Sources

- [configure](../../configure) - target discovery, `--disable-hdf5` filtering, `clean`/`realclean`
- [.github/run_tests.sh](../../.github/run_tests.sh) - broad replay script, downloads, `code_error_workdir`, `pyghl append`, `EXIT` cleanup
- [.github/con2prim-testdata-ref](../../.github/con2prim-testdata-ref) - Con2Prim/Core TestData reference pin
- [.github/et-legacy-testdata-ref](../../.github/et-legacy-testdata-ref) - ET_Legacy/Flux_Source/runner TestData reference pin
- [.github/workflows/github-actions-Ubuntu-clang.yml](../../.github/workflows/github-actions-Ubuntu-clang.yml) - no-HDF5 `c2p_nn_guess` CI leg
- [Unit_Tests/nrpyleakage_main.h](../../Unit_Tests/nrpyleakage_main.h) - key `0`/`1` generation vs. replay
- [Unit_Tests/unit_test_con2prim_tabulated.c](../../Unit_Tests/unit_test_con2prim_tabulated.c) - generation vs. replay modes
- [Unit_Tests/unit_test_c2p_nn_guess.c](../../Unit_Tests/unit_test_c2p_nn_guess.c) - private temporary directory cleanup
- [Unit_Tests/sample_table/generate_simple_table.py](../../Unit_Tests/sample_table/generate_simple_table.py) - checked-in sample-table generator

## See Also

- Depends on: [Build And Configure](../architecture/build-and-configure.md)
- Depends on: [CI Workflows](../architecture/ci-workflows.md)
- See also: [Reference Data And Test Drivers](reference-data-and-test-drivers.md)
- See also: [Contradictions](../contradictions.md#contr-0003)
