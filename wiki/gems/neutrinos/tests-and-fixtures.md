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
- `1`: call `run_unit_test(&eos)` and compare against existing binary fixtures.

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
key `0`; CI and standard unit-test replay use key `1`.

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

No-HDF5 builds set `GHL_DISABLE_HDF5` in `configure`, remove tabulated sources,
and exclude all three NRPyLeakage tests:

- `Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c`
- `Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c`
- `Unit_Tests/unit_test_nrpyleakage_luminosities.c`

For wider build context, see `wiki/build-and-ci.md` and `wiki/test-map.md`.

## Test Behavior Map

`Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c` covers optically thin
gas source replay. Its RHS calls
`NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms`, divides
`R_source` and `Q_source` by `rho`, and evolves `Y_e` and `eps` through RK4.
Each RK4 substep recomputes temperature through `ghl_tabulated_compute_T_from_eps`.
Generation writes the time series; key `1` replays and checks it against
`nrpyleakage_optically_thin_gas_{unperturbed,perturbed}.bin`.

`Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c` covers opacity and
optical-depth behavior. It computes interior/exterior opacities with
`NRPyLeakage_compute_neutrino_opacities`, fills grid arrays, then iterates
optical-depth updates through `NRPyLeakage_optical_depths_PathOfLeastResistance`
until convergence or iteration limit. Key `1` reads and validates
`nrpyleakage_constant_density_sphere_{unperturbed,perturbed}.bin`.

`Unit_Tests/unit_test_nrpyleakage_luminosities.c` has two coverage layers. First,
it checks selected `NRPyLeakage_Fermi_Dirac_integrals` branches for large and
small `z`. Then it replays random metric, primitive, optical-depth, and
luminosity cases through `NRPyLeakage_compute_neutrino_luminosities`, validating
`nue`, `anue`, and `nux` values against
`nrpyleakage_luminosities_{unperturbed,perturbed}.bin`.

`Unit_Tests/unit_test_code_error.c` belongs here only for direct invalid
Fermi-Dirac key coverage: it calls `NRPyLeakage_Fermi_Dirac_integrals(-1, ...)`
for both `z < 1e-3` and `z > 1e-3` cases and maps those keys to
`ghl_error_invalid_fermi_dirac_integral_key`.

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
