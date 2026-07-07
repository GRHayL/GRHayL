# Neutrinos Gem

## Purpose

Neutrinos provides NRPyLeakage support for neutrino opacities, luminosities, optical depths, and GRMHD source terms.

## Read First

- `GRHayL/include/ghl_radiation.h`
- `GRHayL/include/ghl_nrpyleakage.h`
- `GRHayL/include/ghl_nrpyeos_tabulated.h`
- `GRHayL/Neutrinos/NRPyLeakage/`

## Public Headers

- `GRHayL/include/ghl_radiation.h`
- `GRHayL/include/ghl_nrpyleakage.h`
- `GRHayL/include/ghl.h`

Key public surface:
- `ghl_neutrino_luminosities`
- `ghl_neutrino_opacities`
- `ghl_neutrino_optical_depths`
- `NRPyLeakage_Fermi_Dirac_integrals`
- `NRPyLeakage_compute_neutrino_opacities`
- `NRPyLeakage_compute_neutrino_luminosities`
- `NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms`
- `NRPyLeakage_optical_depths_PathOfLeastResistance`

## Implementation Paths

- `GRHayL/Neutrinos/NRPyLeakage/`
- `GRHayL/Neutrinos/make.code.defn`
- `GRHayL/Neutrinos/NRPyLeakage/make.code.defn`

## Test Paths

- `Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c`
- `Unit_Tests/unit_test_nrpyleakage_luminosities.c`
- `Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c`
- `Unit_Tests/nrpyleakage_main.h`

## Key Contracts

- Leakage routines depend on density, `Y_e`, temperature, EOS table access, and optical-depth/opacities structs.
- `ghl_neutrino_optical_depths` is an alias of `ghl_neutrino_opacities`; shape changes affect both concepts.
- Constants and unit conversions are centralized in `GRHayL/include/ghl_nrpyleakage.h`.

## Common Edit Routes

- Change opacity or luminosity formula: update implementation, constants if needed, and targeted leakage unit tests.
- Change radiation structs: update `ghl_radiation.h`, all NRPyLeakage callers, tests, and downstream includes.
- Change EOS dependencies: coordinate with tabulated EOS tests and Con2Prim/primitive variable assumptions.

## Drift Risks

- EOS table variable changes can break leakage even when hydrodynamic tests pass.
- Unit conversion constants are shared by all leakage routines.
- GRHayLib includes radiation and NRPyLeakage headers; public API changes require downstream coordination.

## Do Not Duplicate

Keep formula details and constants in source and public headers. This page gives routing and invariants.
