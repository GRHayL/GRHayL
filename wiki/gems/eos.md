# EOS Gem

## Purpose

EOS provides hybrid, simple, and tabulated equation-of-state support through initialized parameter structs and function pointers.

## Read First

- `docs/raw/EOS.dox`
- `docs/raw/GRHayL_Core.dox`
- `GRHayL/include/ghl.h`
- `wiki/physics/variables-and-conventions.md`

## Public Headers

- `GRHayL/include/ghl.h`
- `GRHayL/include/ghl_eos_functions.h`
- `GRHayL/include/ghl_eos_functions_declaration.h`
- `GRHayL/include/ghl_nrpyeos_hybrid.h`
- `GRHayL/include/ghl_nrpyeos_tabulated.h`

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
