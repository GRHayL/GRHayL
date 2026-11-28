# Neutrinos Tests And Reference Data

> Shared `nrpyleakage_main.h` CLI, binary reference pairs, SLy4 setup, no-HDF5 filtering, per-test coverage, replay tolerances, and CI routes. · Status: confirmed
> Up: [Neutrinos](index.md)

## Summary

NRPyLeakage verification combines table-free analytic/error checks (`unit_test_nrpyleakage_physics.c`, `unit_test_nrpyleakage_classifier_fallback.c`, and `unit_test_code_error.c`) with HDF5 table-backed binary replays (`unit_test_nrpyleakage_optically_thin_gas.c`, `unit_test_nrpyleakage_constant_density_sphere.c`, and `unit_test_nrpyleakage_luminosities.c`) driven through `Unit_Tests/nrpyleakage_main.h`. External BNS_NURATES/DD2 qualification and runtime-cost measurements are defined in [Nucleon Blocking And EOS Conventions](nucleon-blocking-and-eos-conventions.md).

## Detail

### Shared Driver And Binary Reference Pairs

The table-backed NRPyLeakage tests include `Unit_Tests/nrpyleakage_main.h`, which provides the CLI `./test/unit_test_<name> <EOS table path> <test key>` (`0` calls `generate_test_data(&eos)` to write local binaries; `1` calls `run_unit_test(&eos)` to compare against existing binaries; any other key aborts). The driver zero-initializes `eos`, calls `ghl_initialize_tabulated_eos_functions_and_params` with StellarCollapse backend and sound-speed cleaning disabled, sets atmosphere `(rho_b, Y_e, T)`, and calls `ghl_tabulated_free_memory(&eos)` on exit. `unit_test_nrpyleakage_physics.c` and `unit_test_nrpyleakage_classifier_fallback.c` define their own argument-free `main(void)`.

Reference generators live inside the test translation units (no matching `Unit_Tests/data_gen/unit_test_data_nrpyleakage_*.c` files exist) and read or write hard-coded pairs:

- `nrpyleakage_optically_thin_gas_{unperturbed,perturbed}.bin`
- `nrpyleakage_constant_density_sphere_{unperturbed,perturbed}.bin`
- `nrpyleakage_luminosities_{unperturbed,perturbed}.bin`

`.github/run_tests.sh` fetches Neutrinos binaries using `.github/et-legacy-testdata-ref`, whereas compiler workflows fetch the Neutrinos directory from `GRHayL/TestData` `main`. Locally regenerated key-`0` binaries detect implementation drift against installed formulas rather than independently validating physics.

Evidence (descriptive behavior): `nrpyleakage_main.h::main`; corroboration: `run_tests.sh::nrpyleakage_optically_thin_gas`.

### HDF5 Setup And `--disable-hdf5` Build Filtering

Table-backed tests require an unpacked HDF5 EOS table; CI and `.github/run_tests.sh` download `SLy4_3335_rho391_temp163_ye66.h5.bz2`, decompress it to `SLy4_3335_rho391_temp163_ye66.h5`, fetch the binary pairs, and run the table-free executables alongside key-`1` table-backed replays:

```sh
./test/unit_test_nrpyleakage_physics
./test/unit_test_nrpyleakage_classifier_fallback
./test/unit_test_nrpyleakage_optically_thin_gas SLy4_3335_rho391_temp163_ye66.h5 1
./test/unit_test_nrpyleakage_constant_density_sphere SLy4_3335_rho391_temp163_ye66.h5 1
./test/unit_test_nrpyleakage_luminosities SLy4_3335_rho391_temp163_ye66.h5 1
```

Under `--disable-hdf5`, `configure` sets `GHL_DISABLE_HDF5`, retains every guarded NRPyLeakage implementation `.c` file (so `ghl_error_used_disabled_hdf5` return paths remain linkable), excludes `unit_test_nrpyleakage_optically_thin_gas.c`, `unit_test_nrpyleakage_constant_density_sphere.c`, and `unit_test_nrpyleakage_luminosities.c`, and keeps `unit_test_nrpyleakage_physics.c` and `unit_test_nrpyleakage_classifier_fallback.c` built and runnable.

Evidence (descriptive behavior): `configure::unit_test_nrpyleakage_optically_thin_gas`; corroboration: `unit_test_nrpyleakage_classifier_fallback.c::main`.

### Per-Test Coverage And Assertion Strength

- **`Unit_Tests/unit_test_nrpyleakage_physics.c`:** Runs without HDF5 or binary reference files. Compares nucleon blocking and ordinary beta moments against high-precision references across equal, near-equal (cancellation), degenerate, and asymmetric dilute-proton states; checks every Fukushima $F_{-1/2}$ fit interval and dispatcher seam, stable order-zero `log1p`, roundoff-sized fraction normalization, population bounds, transition normalization, zero-shift recovery, threshold orientation, and spectral Kirchhoff pairing.
  Also verifies analytic identities for `6/c` diffusion scaling, `rho^2` bremsstrahlung scaling, the factor of four for heavy leptons in `Q_source`, public-API thick-limit reconstruction, standalone vs. combined opacity equality, single-species endpoint limits, strongly blocked trace-population subnormal ratios (including the luminosity replay's most strongly blocked antineutrino case), all-direction asymmetric optical-depth stencils, exact binary64 bit patterns (infinities, quiet/signaling NaNs, signed zeros, signed minimum subnormals), injected EOS callbacks, output sanitizer fallbacks, and every `ghl_error_used_disabled_hdf5` path under no-HDF5 builds.
- **`Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c`:** Integrates $dY_e/dt = R_{\rm source}/\rho$ and $d\epsilon/dt = Q_{\rm source}/\rho$ via RK4 using `NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms`, recomputing $T$ via `ghl_tabulated_compute_T_from_eps` at each substep. Propagates every EOS inversion status through `ghl_abort_if_error`, tests failure injection before RHS calls and before output file truncation, and compares `time`, `Y_e`, `eps`, and `T` via `ghl_pert_test_fail` (aborting with evolution time on first mismatch).
- **`Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c`:** Computes interior/exterior opacities via `NRPyLeakage_compute_neutrino_opacities` (deriving exterior temperature margin above the table lower bound), iterates `NRPyLeakage_optical_depths_PathOfLeastResistance`, runs `check_sphere_caller_neighbor_order` on the gather helper with unequal x-face metrics, and compares all six opacity and six depth arrays in a serial loop using `ghl_pert_test_fail_with_tolerance` (`1024 * DBL_EPSILON` relative floor, `0.0` absolute tolerance), aborting with grid coordinates on first mismatch.
- **`Unit_Tests/unit_test_nrpyleakage_luminosities.c`:** Checks selected `NRPyLeakage_Fermi_Dirac_integrals` branches (rejecting nonfinite computed or reference values), then compares `nue`, `anue`, and `nux` from `NRPyLeakage_compute_neutrino_luminosities` using `ghl_pert_test_fail_with_tolerance` (`1024 * DBL_EPSILON` relative floor, `0.0` absolute tolerance), aborting with row index on first mismatch. Generation draws each base state once, checks perturbation distance and bounds, and writes row-matched unperturbed/perturbed pairs.
- **`Unit_Tests/unit_test_code_error.c`:** Keys `2` and `3` check `NRPyLeakage_Fermi_Dirac_integrals(-1, ...)` on both `z < 1e-3` and `z > 1e-3`; keys `86` and `87` pass `ghl_error_nrpyleakage_blocking` and `ghl_error_nrpyleakage_nonfinite_output` directly into `expect_error_code` / `ghl_abort_if_error` (without asserting message strings or invoking leakage routines).

Because sphere and luminosity replays set absolute tolerance to `0.0`, tiny opacities, depths, and luminosities are not masked by the generic `1e-30` floor; when both unperturbed and perturbed reference values are exactly `0.0`, the computed value must be `0.0` (a negative control rejects `DBL_MIN`).

Evidence (descriptive behavior): `unit_test_nrpyleakage_constant_density_sphere.c::check_sphere_caller_neighbor_order`; corroboration: `unit_test_nrpyleakage_luminosities.c::luminosity_pert_test_fail`.

### CI Workflow Selection

Every compiler/OS workflow under `.github/workflows/` (Ubuntu GCC, Ubuntu Clang, Ubuntu Intel, macOS GCC, macOS Clang) defines a `neutrinos` job executing the table-free tests and key-`1` replays directly rather than calling `.github/run_tests.sh`. Ubuntu GCC, Ubuntu Clang, and macOS GCC `neutrinos` jobs invoke `.github/actions/code-coverage` (while Intel and macOS Clang comment out the step); for compiler-specific `code-coverage/action.yml` collection behavior and caveats, see [CI Workflows](../../architecture/ci-workflows.md). `.github/run_tests.sh` tracks downloaded/unpacked files in its `EXIT` trap and skips downloads when regular files already exist (see [Runner And Build Outputs](../../tests/runner-and-build-outputs.md)).

Evidence (CI behavior): `github-actions-Ubuntu-gcc.yml::neutrinos`; corroboration: `et-legacy-testdata-ref`.

## Sources

- [Unit_Tests/nrpyleakage_main.h](../../../Unit_Tests/nrpyleakage_main.h) - `main`, key `0`/`1` dispatch, StellarCollapse setup
- [Unit_Tests/unit_test_nrpyleakage_physics.c](../../../Unit_Tests/unit_test_nrpyleakage_physics.c) - table-free blocking, analytic, subnormal, and no-HDF5 tests
- [Unit_Tests/unit_test_nrpyleakage_classifier_fallback.c](../../../Unit_Tests/unit_test_nrpyleakage_classifier_fallback.c) - fallback classifier checks
- [Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c](../../../Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c) - optically thin RK4 replay and lookup-failure injection
- [Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c](../../../Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c) - `check_sphere_caller_neighbor_order`, sphere replay
- [Unit_Tests/unit_test_nrpyleakage_luminosities.c](../../../Unit_Tests/unit_test_nrpyleakage_luminosities.c) - `luminosity_pert_test_fail`, luminosity replay
- [Unit_Tests/unit_test_code_error.c](../../../Unit_Tests/unit_test_code_error.c) - Fermi-Dirac invalid keys `2`/`3`, blocking/nonfinite keys `86`/`87`
- [GRHayL/GRHayL_Core/initialize_eos.c](../../../GRHayL/GRHayL_Core/initialize_eos.c) - `ghl_initialize_tabulated_eos_functions_and_params`
- [configure](../../../configure) - `--disable-hdf5` test exclusions
- [.github/run_tests.sh](../../../.github/run_tests.sh) - local aggregate SLy4/reference download and invocation
- [.github/et-legacy-testdata-ref](../../../.github/et-legacy-testdata-ref) - pinned TestData reference for aggregate runner
- [.github/workflows/github-actions-Ubuntu-gcc.yml](../../../.github/workflows/github-actions-Ubuntu-gcc.yml) - `neutrinos` CI matrix job

## See Also

- Depends on: [Nucleon Blocking And EOS Conventions](nucleon-blocking-and-eos-conventions.md)
- Depends on: [Implementation Flow](implementation-flow.md)
- See also: [Runner And Build Outputs](../../tests/runner-and-build-outputs.md)
- See also: [CI Workflows](../../architecture/ci-workflows.md)
- See also: [Unit Test Catalog](../../tests/unit-test-catalog.md)
