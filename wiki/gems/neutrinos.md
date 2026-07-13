# Neutrinos Gem

## Purpose

This page routes Neutrinos KB work. It summarizes where to read first for
NRPyLeakage public data, implementation flow, and tests without replacing
source, headers, tests, CI, or Doxygen source.

## Read Order

1. [API And Data](neutrinos/api-and-data.md) for radiation structs, public
   `NRPyLeakage_*` declarations, constants ownership, errors, and HDF5/EOS
   dependency.
2. [Implementation Flow](neutrinos/implementation-flow.md) for the five
   `GRHayL/Neutrinos/NRPyLeakage/` source files and their writeback paths.
3. [Tests And Fixtures](neutrinos/tests-and-fixtures.md) for unit tests,
   fixture pairs, EOS table setup, and CI downloads.

## Ground Truth

- Source: `GRHayL/Neutrinos/NRPyLeakage/`
- Public radiation structs: `GRHayL/include/ghl_radiation.h`
- Public leakage declarations and constants: `GRHayL/include/ghl_nrpyleakage.h`
- Error codes and EOS parameter types: `GRHayL/include/ghl.h`
- Direct tabulated EOS dependencies: `GRHayL/include/ghl_nrpyeos_tabulated.h`
  and `GRHayL/include/ghl_eos_functions.h`
- Tests: `Unit_Tests/nrpyleakage_main.h` and
  `Unit_Tests/unit_test_nrpyleakage_*.c`
- Fermi-Dirac error coverage: `Unit_Tests/unit_test_code_error.c`
- Build and CI gates: `configure`, `.github/run_tests.sh`, and
  `.github/workflows/`

If this page or a child page conflicts with those files, trust the repo-local
ground truth and update the KB.

## Public Surface

- `ghl_neutrino_luminosities`
- `ghl_neutrino_opacities`
- `ghl_neutrino_optical_depths`
- `NRPyLeakage_Fermi_Dirac_integrals`
- `NRPyLeakage_compute_neutrino_opacities`
- `NRPyLeakage_compute_neutrino_luminosities`
- `NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms`
- `NRPyLeakage_optical_depths_PathOfLeastResistance`

The `NRPyLeakage_*` spelling is the exported family; no `ghl_*` wrapper family
exists for these calls. Radiation container types retain the `ghl_` prefix.

## Contract Summary

- Radiation struct order is `nue`, `anue`, `nux`; opacity/depth has two slots
  per species. Slot meanings are used differently by number and energy source
  formulas but remain unnamed in the public header, so keep `[0]`/`[1]`
  wording in caller contracts.
- EOS-dependent routines expect initialized HDF5-backed tabulated EOS state,
  geometric density, EOS-compatible positive temperature, and valid output
  pointers. They return before writeback on HDF5, EOS, or generated Fermi
  errors.
- Optical-depth update is a `void` six-neighbor stencil call with no validation
  or failure channel. Metric stencil order is minus/center/plus.
- All five implementation files match their manifest and header declarations.
- The aggregate runner invokes all three HDF5 test binaries, and all five
  workflow families configure those commands. Command presence is not a
  historical execution result. When executed, their main fixture comparisons
  discard `ghl_pert_test_fail` results, so a completed run establishes setup,
  execution, and fixture reads but not successful numerical replay.

## Scope Notes

- Neutrinos KB pages live under `wiki/` only; Doxygen source under `docs/**`
  is a separate authority (currently no dedicated Neutrinos page exists there).
- Keep generated formula blocks in `GRHayL/Neutrinos/NRPyLeakage/*.c`.
- Keep constants and unit conversions in `GRHayL/include/ghl_nrpyleakage.h`;
  do not copy them into KB tables.
- Treat HDF5-enabled tabulated EOS access as part of the direct Neutrinos
  dependency surface because public leakage routines call it.
