# Neutrinos API And Data

## Purpose

This page maps Neutrinos public data shapes, callable entry points, constants,
error behavior, and HDF5/EOS dependencies. Ground truth remains the named
repo-local headers, source files, tests, and build scripts.

## Radiation Data

`GRHayL/include/ghl_radiation.h` owns the public radiation structs:

- `ghl_neutrino_luminosities` has scalar `nue`, `anue`, and `nux` fields for
  electron neutrino, electron antineutrino, and heavy-lepton neutrino
  luminosities.
- `ghl_neutrino_opacities` has `nue[2]`, `anue[2]`, and `nux[2]`.
- `ghl_neutrino_optical_depths` is a typedef alias of
  `ghl_neutrino_opacities`, so depth and opacity storage share the same shape.

Use `[0]` and `[1]` as two entries per species unless implementation source
proves stronger semantics for a specific call path. Current public headers name
the arrays but do not publish semantic labels for those two entries.

`ghl_radiation.h` includes `ghl.h` first, defines the radiation structs, then
includes `ghl_nrpyleakage.h`; leakage declarations can therefore use the
radiation types.

## Public NRPyLeakage Calls

`GRHayL/include/ghl_nrpyleakage.h` declares five public routines:

- `NRPyLeakage_Fermi_Dirac_integrals`
- `NRPyLeakage_compute_neutrino_opacities`
- `NRPyLeakage_compute_neutrino_luminosities`
- `NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms`
- `NRPyLeakage_optical_depths_PathOfLeastResistance`

The first four return `ghl_error_codes_t`. The optical-depth path routine
returns `void` and writes `ghl_neutrino_optical_depths` output.

## Constants And Units

`GRHayL/include/ghl_nrpyleakage.h` owns leakage switches, physical constants,
unit conversions, and derived constants. KB pages should route to that header
instead of duplicating a constants table. Formula blocks that consume those
constants live in `GRHayL/Neutrinos/NRPyLeakage/*.c`.

The same header defines `robust_isnan`, `robust_isfinite`, and
`NRPYLEAKAGE_FD_OR_RETURN`. That macro calls
`NRPyLeakage_Fermi_Dirac_integrals` and returns the error immediately when the
helper fails.

## Fermi-Dirac Error Behavior

`GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_Fermi_Dirac_integrals.c` supports
keys `0` through `5` and splits approximations on `z > 1e-3`. It zeroes the
output value before branch dispatch. Unsupported keys return
`ghl_error_invalid_fermi_dirac_integral_key`, defined in
`GRHayL/include/ghl.h`.

`Unit_Tests/unit_test_code_error.c` exercises invalid Fermi-Dirac keys on both
`z` branches and maps them to `ghl_error_invalid_fermi_dirac_integral_key`.

## HDF5 And EOS

Neutrinos source routines use tabulated EOS calls directly:

- `NRPyLeakage_compute_neutrino_opacities`
- `NRPyLeakage_compute_neutrino_luminosities`
- `NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms`

Those three routines call
`ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T` to obtain chemical
potentials and composition. Under `GHL_DISABLE_HDF5`, each returns
`ghl_error_used_disabled_hdf5` before table access.

The optical-depth path routine does not call EOS or HDF5. It consumes local and
neighbor `ghl_neutrino_opacities`/`ghl_neutrino_optical_depths`, metric stencil
arrays, and grid spacing, then writes one optical-depth struct.

Neutrinos tests use additional EOS setup calls that are test-only dependencies
for fixture generation or replay:

- `Unit_Tests/nrpyleakage_main.h` calls
  `ghl_initialize_tabulated_eos_functions_and_params` from `GRHayL/include/ghl.h`.
- `Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c` calls
  `ghl_tabulated_compute_eps_from_T` and `ghl_tabulated_compute_T_from_eps`
  from `GRHayL/include/ghl_eos_functions.h`.

`configure` adds `GHL_DISABLE_HDF5` and excludes the three
`unit_test_nrpyleakage_*.c` tests when HDF5 is disabled. `.github/run_tests.sh`
downloads the SLy4 EOS table and Neutrinos fixture pairs before running the
three NRPyLeakage unit tests with key `1`.
