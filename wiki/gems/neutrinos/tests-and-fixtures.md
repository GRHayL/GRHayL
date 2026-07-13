# Neutrinos Tests And Fixtures

Purpose: route NRPyLeakage test setup, fixture data, and CI commands back to
repo-local ground truth. Source, tests, `configure`, and CI scripts remain the
authority.

## Shared Harness

All three NRPyLeakage unit-test files include `Unit_Tests/nrpyleakage_main.h`.
That helper owns the common CLI:

```sh
./test/unit_test_<nrpyleakage test> <EOS table path> <test key>
```

The accepted keys are:

- `0`: call `generate_test_data(&eos)` and write local binary reference data.
- `1`: call `run_unit_test(&eos)` and attempt comparison against existing
  binary fixtures; see the assertion gap below.

Any other key is rejected by the shared harness. The harness also initializes a
tabulated, stellar-collapse EOS before dispatching either mode. Note: the
harness pre-sets `eos.eos_type`, `eos.table_type`, and
`eos.clean_sound_speed = true`, but the subsequent
`ghl_initialize_tabulated_eos_functions_and_params(...)` call overwrites all
three -- it hard-codes `ghl_eos_table_stellarcollapse` and
`clean_sound_speed = false` (`GRHayL/GRHayL_Core/initialize_eos.c`), so the
effective EOS runs with sound-speed cleaning disabled. The harness then
passes atmosphere values for `rho_b`, `Y_e`, and `T`, and frees table memory
with `ghl_tabulated_free_memory(&eos)` after either mode completes.

Ground truth:

- `Unit_Tests/nrpyleakage_main.h`
- `GRHayL/include/ghl_nrpyleakage.h`
- `GRHayL/include/ghl_radiation.h`

## Fixture Pairs

Fixture names are hard-coded by each test and are downloaded by
`.github/run_tests.sh` for normal CI replay:

- `nrpyleakage_optically_thin_gas_unperturbed.bin`
- `nrpyleakage_optically_thin_gas_perturbed.bin`
- `nrpyleakage_constant_density_sphere_unperturbed.bin`
- `nrpyleakage_constant_density_sphere_perturbed.bin`
- `nrpyleakage_luminosities_unperturbed.bin`
- `nrpyleakage_luminosities_perturbed.bin`

The same set can be written locally by running the corresponding executable with
key `0`; CI and standard unit-test replay use key `1`. These producers live in
the three test translation units themselves; there are no matching
`Unit_Tests/data_gen/unit_test_data_nrpyleakage_*.c` targets.

## HDF5 And EOS Setup

NRPyLeakage tests require an HDF5-backed tabulated EOS path. The repo-local CI
route downloads `SLy4_3335_rho391_temp163_ye66.h5.bz2`, unpacks it to
`SLy4_3335_rho391_temp163_ye66.h5`, downloads the six Neutrinos fixture files,
then runs:

```sh
./test/unit_test_nrpyleakage_optically_thin_gas SLy4_3335_rho391_temp163_ye66.h5 1
./test/unit_test_nrpyleakage_constant_density_sphere SLy4_3335_rho391_temp163_ye66.h5 1
./test/unit_test_nrpyleakage_luminosities SLy4_3335_rho391_temp163_ye66.h5 1
```

Use `.github/run_tests.sh` as the primary command route for fixture and table
downloads. Do not claim current remote fixture availability from this page; the
repo-local script is the cited route.

No-HDF5 builds set `GHL_DISABLE_HDF5` and apply `configure`'s exact
implementation-source filter. That filter removes many tabulated EOS/Flux_Source
sources but retains selected Con2Prim Tabulated/NN helpers and all NRPyLeakage
implementation files. It also excludes all three NRPyLeakage tests:

- `Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c`
- `Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c`
- `Unit_Tests/unit_test_nrpyleakage_luminosities.c`

NRPyLeakage implementation sources remain compiled so their early
`ghl_error_used_disabled_hdf5` paths exist, but exclusion of all three tests
means those paths have no direct no-HDF5 test. `unit_test_code_error` keys `2`
and `3` cover only invalid Fermi keys and remain available without HDF5.

For wider build context, see `wiki/build-and-ci.md` and `wiki/test-map.md`.

## Test Behavior Map

`Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c` covers optically thin
gas source replay. Its RHS calls
`NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms`, divides
`R_source` and `Q_source` by `rho`, and evolves `Y_e` and `eps` through RK4.
Each RK4 substep recomputes temperature through `ghl_tabulated_compute_T_from_eps`.
Generation writes the time series; key `1` replays and calls comparison helpers
against
`nrpyleakage_optically_thin_gas_{unperturbed,perturbed}.bin`.

`Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c` covers opacity and
optical-depth behavior. It computes interior/exterior opacities with
`NRPyLeakage_compute_neutrino_opacities`, fills grid arrays, then iterates
optical-depth updates through `NRPyLeakage_optical_depths_PathOfLeastResistance`
until convergence or iteration limit. Key `1` reads and calls comparison
helpers on
`nrpyleakage_constant_density_sphere_{unperturbed,perturbed}.bin`.

`Unit_Tests/unit_test_nrpyleakage_luminosities.c` has two coverage layers. First,
it checks selected `NRPyLeakage_Fermi_Dirac_integrals` branches for large and
small `z`. Then it replays random metric, primitive, optical-depth, and
luminosity cases through `NRPyLeakage_compute_neutrino_luminosities`, then calls
comparison helpers for `nue`, `anue`, and `nux` values against
`nrpyleakage_luminosities_{unperturbed,perturbed}.bin`.

`Unit_Tests/unit_test_code_error.c` belongs here only for direct invalid
Fermi-Dirac key coverage: it calls `NRPyLeakage_Fermi_Dirac_integrals(-1, ...)`
for both `z < 1e-3` and `z > 1e-3` cases and maps those keys to
`ghl_error_invalid_fermi_dirac_integral_key`.

### Effective Assertion Gap

All three NRPyLeakage replay tests call `ghl_pert_test_fail`, which returns a
boolean, but discard that return value instead of branching to `ghl_error`.
Consequences:

- Optically thin replay reads and recomputes time, `Y_e`, `eps`, and `T`, but a
  numerical mismatch cannot fail the executable.
- Constant-density-sphere replay reads and recomputes all six opacity and six
  optical-depth arrays, but a numerical mismatch cannot fail the executable.
- Luminosity replay reads and recomputes `nue`, `anue`, and `nux`, but a
  numerical mismatch cannot fail the executable.

File-open/read failures, EOS/leakage errors passed to `ghl_abort_if_error`, and
the luminosity test's explicit valid Fermi checks can still fail. Classify all
three as compile + execute + fixture-read evidence with ineffective main
numerical replay assertions, not validated numerical replay.

Further bounded gaps:

- Luminosity explicitly checks only high-branch Fermi keys `0`, `1`, `2` and
  low-branch keys `0`, `1`; valid keys `3` through `5` are indirect only.
- Constant-density-sphere passes every optical-depth neighbor pair in
  plus-then-minus order, opposite the function signature. Flat metric stencils
  and minimum selection hide this reversal, so unequal-direction mapping is
  untested.
- Constant-density-sphere also passes perturbed fixture values as the
  `trusted` argument and unperturbed values as the `perturbed` argument to
  `ghl_pert_test_fail`, opposite the helper contract. Discarded return values
  already make these calls ineffective; fixing only the discard would leave
  this argument-order fault.
- Combined source-term opacity outputs are computed but never compared in the
  optically thin test.
- No Neutrinos test deliberately exercises EOS-error propagation,
  disabled-HDF5 returns, non-finite final luminosity/source outputs, or the
  `EnsureFinite`/`robust_isfinite` fallback paths.

## CI Routes

`.github/run_tests.sh` is the compact end-to-end route for Neutrinos fixture
downloads and test commands. The concrete workflow example is
`.github/workflows/github-actions-Ubuntu-gcc.yml`: its `neutrinos` job matrix
lists `nrpyleakage_optically_thin_gas`,
`nrpyleakage_constant_density_sphere`, and `nrpyleakage_luminosities`, then
downloads matching `_unperturbed.bin` and `_perturbed.bin` files plus the SLy4
table before running each test with key `1`.

Matching Neutrinos jobs exist in other workflow files. Before expanding line
citations beyond the Ubuntu GCC example, verify exact locations with:

```sh
rg -n "neutrinos|nrpyleakage" .github/workflows
```

Exact workflow set is five files: Ubuntu GCC, Ubuntu Clang, Ubuntu Intel,
macOS GCC, and macOS Clang. Each uses two OS versions and all three test names.
Coverage upload is active for GCC and Ubuntu Clang Neutrinos jobs; it is
commented out for Intel and macOS Clang. This does not repair the discarded
comparison results.

The aggregate runner deletes working-directory `*.bin`, `*.h5`, and `*.bz2`
files at the end. Individual workflow jobs rely on fresh workspaces and do not
include an explicit Neutrinos cleanup step.

## Repo-Local References

- `Unit_Tests/nrpyleakage_main.h`
- `Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c`
- `Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c`
- `Unit_Tests/unit_test_nrpyleakage_luminosities.c`
- `Unit_Tests/unit_test_code_error.c`
- `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c`
- `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities.c`
- `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_luminosities.c`
- `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_Fermi_Dirac_integrals.c`
- `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_optical_depths_PathOfLeastResistance.c`
- `GRHayL/include/ghl_nrpyleakage.h`
- `GRHayL/include/ghl_radiation.h`
- `configure`
- `.github/run_tests.sh`
- `.github/workflows/github-actions-Ubuntu-gcc.yml`
- `wiki/build-and-ci.md`
- `wiki/test-map.md`
