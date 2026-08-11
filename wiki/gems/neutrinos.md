# Neutrinos Gem

## Purpose

This page routes Neutrinos KB work. It summarizes where to read first for
NRPyLeakage physics, generated-code provenance, public data, implementation
flow, portability, and tests without replacing source, headers, tests, CI, or
Doxygen source.

## Read Order

1. [Physics And EOS Contract](neutrinos/physics-and-eos-contract.md) for
   species, chemical potentials, free-nucleon fractions, units, and the
   number and energy source-term conventions.
2. [Generator Provenance](neutrinos/generator-provenance.md) for the original
   Python notebooks, symbolic common-subexpression temporaries, and the
   boundary between ancestral generators and current checked-in C.
3. [API And Data](neutrinos/api-and-data.md) for radiation structs, public
   `NRPyLeakage_*` declarations, constants ownership, errors, and HDF5/EOS
   dependency.
4. [CompOSE EOS Adapter How-To](neutrinos/compose-eos-adapter-how-to.md) for
   adapting CompOSE state, composition, and chemical-potential outputs to the
   NRPyLeakage EOS callback and validating that boundary.
5. [Implementation Flow](neutrinos/implementation-flow.md) for the five
   `GRHayL/Neutrinos/NRPyLeakage/` source files, their writeback paths, and
   the minimal direct-compilation boundary.
6. [Tests And Fixtures](neutrinos/tests-and-fixtures.md) for unit tests,
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
- Ancestral derivation and generator evidence: the
  [Tabulated_EOS_IllinoisGRMHD repository](https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD),
  with durable roles and conventions summarized in
  [Generator Provenance](neutrinos/generator-provenance.md) and
  [Physics And EOS Contract](neutrinos/physics-and-eos-contract.md)

If this page, a child page, or external evidence conflicts with current
repo-local source, headers, manifests, or tests, trust the repo-local ground
truth and update the KB.

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
- The external notebooks are provenance, not the authority for current GRHayL
  signatures or behavior. Preserve the durable derivation and interface facts
  in child pages so routine KB use does not depend on that repository remaining
  available.
- Keep constants and unit conversions in `GRHayL/include/ghl_nrpyleakage.h`;
  do not copy them into KB tables.
- Treat HDF5-enabled tabulated EOS access as part of the direct dependency
  surface for the opacity, combined source/opacity, and luminosity routines;
  the Fermi helper and optical-depth update do not call the EOS.

## Ground Truth References

- [Original Tabulated_EOS_IllinoisGRMHD repository](https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD)
