# Con2Prim Tests And Reference Data

> Hybrid, tabulated, neural-network, limit, conversion, ET_Legacy, and debug Con2Prim unit tests, sanitizer requirements, reference data, and targeted commands. · Status: confirmed
> Up: [Con2Prim](index.md)

## Summary

Con2Prim tests span shared-data hybrid/limit/conversion replay (`unit_test_con2prim_multi_method_hybrid.c` and the limit/conversion drivers, replaying `unit_test_data_con2prim_multi_method_hybrid.c` output), hard-coded failure/round-trip checks (`unit_test_hybrid_failure.c`), HDF5-backed tabulated replay (`unit_test_con2prim_tabulated.c` with `SLy4` and `pyghl append`), self-contained NN checks (`unit_test_c2p_nn_guess.c`), ET_Legacy replay, and manual-only `unit_test_con2prim_debug.c`. Newman final-iteration memory safety requires ASan/UBSan.

## Detail

### Hybrid, Limit, Conversion, And ET_Legacy Tests

| Test binary | Coverage scope & key invariants | Reference data & runner mode |
| --- | --- | --- |
| `unit_test_con2prim_multi_method_hybrid.c` | Hybrid selected-method recovery (`Font1D`, `Palenzuela1D`, `Noble1D_entropy`, `Palenzuela1D_entropy`, `Noble1D`, `Noble2D`) passing `cons_undens` to guess and selector; Noble negative-pressure cases near `P ~ 1e-30` must lie at roundoff scale relative to enthalpy density. (Ev: `unit_test_con2prim_multi_method_hybrid.c`) | `metric_Bfield_initial_data.bin`, `con2prim_multi_method_hybrid_*.bin` (`data_gen/unit_test_data_con2prim_multi_method_hybrid.c`); default runner & `c2p-routines`. |
| `unit_test_hybrid_failure.c` | Hard-coded `Noble2D` max-iter/singular failures, finite-difference Noble derivatives across unequal-Gamma pieces, independent Newton budgets, direct `ghl_brent` checks, Font1D curved/multi-piece/cyclic-EOS round trips, `utilde` monotone rescaling, and `Noble1D_entropy2` round trips and incompatible-EOS rejection. (Ev: `unit_test_hybrid_failure.c`) | Hard-coded (no external data); default runner & `c2p-failure`. |
| `unit_test_apply_conservative_limits.c` | Poisoned diagnostics init, low-field and high-`psi6` momentum rescaling, and synthetic energies around `tau_atm`. (Ev: `unit_test_apply_conservative_limits.c`) | Shared `metric_Bfield_initial_data.bin` + `apply_conservative_limits_*.bin`; default runner & `c2p-routines`. |
| `unit_test_enforce_primitive_limits_and_compute_u0.c` | Wrapper floors/ceilings, entropy recompute, and `u0`, plus direct `ghl_limit_v_and_compute_u0` cap/non-limit cases. (Ev: `unit_test_enforce_primitive_limits_and_compute_u0.c`) | Shared `metric_Bfield_initial_data.bin` + `enforce_primitive_limits_and_compute_u0_*.bin`; default runner & `c2p-routines`. |
| `unit_test_compute_conservs_and_Tmunu.c` & `unit_test_ET_Legacy_{conservs,primitives}.c` | Forward conservative/`T_{\mu\nu}` checks and ET_Legacy compatibility replay. (Ev: `unit_test_ET_Legacy_primitives.c`) | Shared `metric_Bfield_initial_data.bin` + `compute_conservs_and_Tmunu_*.bin`; `ET_Legacy_*.bin`; default runner & `c2p-routines`/`ET-Legacy`. |

### Tabulated, NN Guess, Debug, And Sanitizer Requirements

- **`unit_test_con2prim_tabulated.c` (HDF5-only):** Replays `Palenzuela1D`, `Newman1D_entropy`, and `Newman1D`/`Palenzuela1D_entropy`/`Noble2D` with backups over `rho_vs_T` and `Pmag_vs_Wm1` grids using `SLy4_3335_rho391_temp163_ye66.h5` (prepared via `pyghl append` in runner/CI before key `1`; key `0` generates local data). Replays NN disabled and enabled, directly reproducing every successful NN main retry. Also covers speed-limited Noble enthalpy inversion, fixed-temperature Palenzuela (without exercising failed-Brent `T_min` reset), tabulated Noble callback error paths, and Newman cubic/nonfinite errors.
- **Newman history sanitizer requirement (`Unit_Tests/con2prim_newman_history.h`):** Drives tabulated Newman energy and entropy helpers to iteration exhaustion (`600` EOS calls = `maxsteps = 300` x 2 with `T_min` retry) without Aitken resets. Ordinary exit/call-count assertions do not verify memory safety; detecting an out-of-bounds history write requires building/running with ASan/UBSan.
- **`unit_test_c2p_nn_guess.c`:** Verifies `GHL_NN_C2P_API_VERSION` layout, overflow-safe `x` fallback, atmosphere fallback, root/embedded/legacy HDF5 loading, malformed datasets, and failed-load preservation in a private temporary directory (HDF5 setup compiled out under `GHL_DISABLE_HDF5`; tested in Ubuntu-Clang `c2p-failure` no-HDF5 leg).
- **`unit_test_con2prim_debug.c`:** Manual-only single-case tabulated `Newman1D` stderr dump taking a CLI table path; excluded from `run_tests.sh` and `--disable-hdf5`.

Evidence (descriptive behavior): `unit_test_con2prim_tabulated.c`; corroboration: `con2prim_newman_history.h` and `unit_test_c2p_nn_guess.c`.

### Targeted Reruns And Companion TestData Coordination

In a disposable checkout, after `./configure -r && make tests` and exporting `LD_LIBRARY_PATH="$(pwd)/build/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"`, run only the binaries matching the edit from the repository root as `./test/<name>`. `unit_test_hybrid_failure` and `unit_test_c2p_nn_guess` need no external data. `unit_test_apply_conservative_limits`, `unit_test_con2prim_multi_method_hybrid`, `unit_test_enforce_primitive_limits_and_compute_u0`, `unit_test_compute_conservs_and_Tmunu`, and `unit_test_ET_Legacy_{conservs,primitives}` read the `.bin` files listed in the `Reference data & runner mode` column from the working directory, so first stage the `GRHayL/TestData` files as [Runner And Build Outputs](../../tests/runner-and-build-outputs.md#githubrun_testssh-reference-selectors-and-matrix-differences) describes.

`./test/unit_test_con2prim_tabulated SLy4_3335_rho391_temp163_ye66.h5 1` also needs the staged `con2prim_tabulated_*` files and a table prepared with `pyghl append SLy4_3335_rho391_temp163_ye66.h5`, because key `1` loads the `grhayl_nn_c2p` model group from that table.

Generated hybrid data are implementation-derived regression bars (the forward round trips in `unit_test_hybrid_failure.c` supply independent closure checks). Keep producer and reader synchronized for entropy settings and never alter signed perturbations or the field-aligned momentum edge merely to widen tolerances. When updating companion `GRHayL/TestData` files, publish an immutable revision before updating `.github/con2prim-testdata-ref`.

## Sources

- [Unit_Tests/unit_test_con2prim_multi_method_hybrid.c](../../../Unit_Tests/unit_test_con2prim_multi_method_hybrid.c) - hybrid selected-method replay
- [Unit_Tests/data_gen/unit_test_data_con2prim_multi_method_hybrid.c](../../../Unit_Tests/data_gen/unit_test_data_con2prim_multi_method_hybrid.c) - shared hybrid/limit/Core generator
- [Unit_Tests/unit_test_hybrid_failure.c](../../../Unit_Tests/unit_test_hybrid_failure.c) - hard-coded hybrid failure, Brent, Font1D, and `Noble1D_entropy2` checks
- [Unit_Tests/unit_test_con2prim_tabulated.c](../../../Unit_Tests/unit_test_con2prim_tabulated.c) - tabulated multi-method and NN replay
- [Unit_Tests/con2prim_newman_history.h](../../../Unit_Tests/con2prim_newman_history.h) - final-iteration Newman callback helper
- [Unit_Tests/unit_test_c2p_nn_guess.c](../../../Unit_Tests/unit_test_c2p_nn_guess.c) - NN inference and HDF5 loader checks
- [Unit_Tests/unit_test_con2prim_debug.c](../../../Unit_Tests/unit_test_con2prim_debug.c) - manual tabulated debug driver
- [Unit_Tests/unit_test_apply_conservative_limits.c](../../../Unit_Tests/unit_test_apply_conservative_limits.c) - conservative limit tests
- [Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c](../../../Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c) - primitive limit and `u0` tests
- [Unit_Tests/unit_test_ET_Legacy_primitives.c](../../../Unit_Tests/unit_test_ET_Legacy_primitives.c) - ET_Legacy C2P replay
- [.github/run_tests.sh](../../../.github/run_tests.sh) - `con2prim/` downloads and `pyghl append`

## See Also

- Depends on: [Reference Data And Test Drivers](../../tests/reference-data-and-test-drivers.md)
- Depends on: [Data Generators And Helpers](../../tests/data-generators-and-helpers.md)
- See also: [Solver Matrix](solver-matrix.md)
- See also: [Neural-Network Primitive Guess](neural-network-primitive-guess.md)
