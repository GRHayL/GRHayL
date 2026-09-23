# EOS Gem

## Purpose

EOS provides hybrid, simple, and tabulated equation-of-state support through initialized parameter structs and function pointers.

## Read First

- [EOS initialization and dispatch](eos/initialization-and-dispatch.md)
- [Hybrid piecewise-polytrope EOS](eos/hybrid-piecewise-polytrope.md)
- [Tabulated EOS table contract](eos/tabulated-table-contract.md)
- [Stellarcollapse table adapter](eos/stellarcollapse-table-adapter.md)
- [Tabulated EOS interpolation and bounds](eos/tabulated-interpolation-and-bounds.md)
- [Tabulated interpolator catalog](eos/tabulated-interpolator-catalog.md)
- [CompOSE-to-NRPyLeakage adapter how-to](neutrinos/compose-eos-adapter-how-to.md)
- [EOS tests and fixtures](eos/tests-and-fixtures.md)
- `docs/raw/EOS.dox`
- `docs/raw/GRHayL_Core.dox`
- `GRHayL/include/ghl.h`
- `wiki/physics/variables-and-conventions.md`

## Public Headers

- `GRHayL/include/ghl.h`
- `GRHayL/include/ghl_eos_functions.h`
- `GRHayL/include/ghl_nrpyeos_hybrid.h`
- `GRHayL/include/ghl_nrpyeos_tabulated.h`

`GRHayL/include/ghl_eos_functions_declaration.h` is an internal,
non-installed storage-definition companion.

Key public surface:
- `ghl_initialize_simple_eos_functions_and_params`
- `ghl_initialize_hybrid_eos_functions_and_params`
- `ghl_initialize_tabulated_eos_functions_and_params`
- `ghl_compute_h_and_cs2`
- Hybrid `ghl_hybrid_*` function pointers
- Tabulated `ghl_tabulated_*` function pointers

## Implementation Paths

- Hybrid EOS: `GRHayL/EOS/Hybrid/`
- Tabulated EOS: `GRHayL/EOS/Tabulated/`
- Tabulated interpolators: `GRHayL/EOS/Tabulated/interpolators/`
- Stellar-collapse table adapter: `GRHayL/EOS/Tabulated/stellarcollapse/`
- Core initialization: `GRHayL/GRHayL_Core/initialize_eos.c`

Focused routes:
- Simple/ideal-fluid setup and EOS function-pointer dispatch:
  [initialization and dispatch](eos/initialization-and-dispatch.md).
- Hybrid/piecewise-polytrope fields, cold pressure, cold energy, entropy,
  `K_ppoly`, and `eps_integ_const`:
  [hybrid piecewise-polytrope EOS](eos/hybrid-piecewise-polytrope.md).
- Tabulated HDF5 table lifecycle, sample table roles, energy shift, and memory:
  [tabulated table contract](eos/tabulated-table-contract.md).
- Stellar-collapse table-type dispatch, HDF5 helper boundary, dataset
  ownership, and conversion into `ghl_eos_parameters`:
  [stellarcollapse table adapter](eos/stellarcollapse-table-adapter.md).
- Tabulated bounds, interpolators, temperature recovery, and beta equilibrium:
  [tabulated interpolation and bounds](eos/tabulated-interpolation-and-bounds.md).
- Built tabulated interpolator wrappers by input tuple, beta-equilibrium/rho
  maps, enthalpy tabulation, and sound-speed cleaning:
  [tabulated interpolator catalog](eos/tabulated-interpolator-catalog.md).
- CompOSE-to-NRPyLeakage EOS mapping, native-adapter boundary, integration
  routes, and validation: [CompOSE EOS adapter how-to](neutrinos/compose-eos-adapter-how-to.md).
- EOS unit tests, helper-only files, sample tables, generated `simple_table.h5`,
  and no-HDF5 error behavior: [tests and fixtures](eos/tests-and-fixtures.md).

## Test Paths

- `Unit_Tests/unit_test_piecewise_polytrope.c`
- `Unit_Tests/unit_test_tabulated_eos.c`
- `Unit_Tests/test_compute_h_and_cs2.c`
- `Unit_Tests/tabulated_eos_unit_test_helpers.c`
- `Unit_Tests/sample_table/`

## Key Contracts

- `ghl_eos_parameters` stores EOS type, table type, atmosphere values, bounds, hybrid pieces, and tabulated table state.
- Initialize EOS through core initialization helpers so function pointers match EOS type.
- Hybrid EOS splits cold and thermal behavior; tabulated EOS interpolates table quantities and enforces table bounds.
- Tabulated EOS runtime table support depends on HDF5-enabled builds.
- `rho`, pressure, internal energy, entropy, temperature, and `Y_e` bounds must stay consistent with primitive and conservative limiters.

Current support boundaries:

- A live tabulated EOS must be explicitly cleaned before reinitialization or an
  EOS-family switch; live-object replacement is unsupported.
- Tabulated `tau_atm` uses minimum density/energy rather than its atmosphere
  density/energy pair.
- Table validation establishes the structural and indexing domain used by the
  interpolators; it is not a complete physical or malformed-HDF5 audit.

Route exact evidence through initialization, hybrid, table, and interpolator
leaves above; do not broaden these statements into guarantees outside the
documented lifecycle and input domains.

## Common Edit Routes

- Add hybrid calculation: update `GRHayL/EOS/Hybrid/`, public pointer/declaration if exported, initialization assignment, tests, and `docs/raw/EOS.dox`.
- Add tabulated interpolator: update `GRHayL/EOS/Tabulated/interpolators/`, `ghl_nrpyeos_tabulated.h`, function pointer files if public, table tests, and build lists.
- Change EOS initialization: update `GRHayL/GRHayL_Core/initialize_eos.c`, core docs, and tests that initialize EOS.
- Change table layout support: update table adapter, HDF5 helpers, sample-table tests, and no-HDF5 guard behavior.

## Drift Risks

- Function pointer declarations, definitions, and initialization assignments can drift apart.
- Con2Prim, Flux_Source, and Neutrinos depend on EOS variables beyond pressure, especially entropy, temperature, and `Y_e`.
- Downstream GRHayLib includes EOS public headers and compiles EOS subdirectories; source-tree shape changes require coordination.

## Do Not Duplicate

Doxygen and headers are authority for EOS APIs and table details. Keep long table-variable lists and interpolation details in source docs.
