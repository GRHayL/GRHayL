# Flux Source Gem

## Purpose

Flux_Source computes characteristic speeds, HLLE fluxes, and hydrodynamic source terms used by the conservative evolution equations.

## Read First

- [Characteristic speeds contract](flux-source/characteristic-speeds-contract.md)
- [HLLE flux variant matrix](flux-source/hlle-flux-variant-matrix.md)
- [Source-term contract](flux-source/source-terms-contract.md)
- [Generated NRPy boundary](flux-source/generated-nrpy-boundary.md)
- `docs/raw/Flux_Source.dox`
- `docs/raw/derivation.md`
- `docs/raw/Reconstruction.dox`
- `GRHayL/include/ghl_flux_source.h`
- `wiki/physics/evolution-equation-map.md`

## Public Headers

- `GRHayL/include/ghl_flux_source.h`
- `GRHayL/include/ghl.h`

Key public surface:
- `ghl_calculate_source_terms`
- `ghl_calculate_characteristic_speed_dirn0`
- `ghl_calculate_characteristic_speed_dirn1`
- `ghl_calculate_characteristic_speed_dirn2`
- `ghl_calculate_HLLE_fluxes_dirn*_hybrid`
- `ghl_calculate_HLLE_fluxes_dirn*_hybrid_entropy`
- `ghl_calculate_HLLE_fluxes_dirn*_tabulated`
- `ghl_calculate_HLLE_fluxes_dirn*_tabulated_entropy`

## Implementation Paths

- Shared speeds: `GRHayL/Flux_Source/`; route through
  [characteristic speeds contract](flux-source/characteristic-speeds-contract.md).
- Source terms: `GRHayL/Flux_Source/ghl_calculate_source_terms.c`; route
  through [source-term contract](flux-source/source-terms-contract.md).
- Hybrid, hybrid entropy, tabulated, and tabulated entropy fluxes: route
  through [HLLE flux variant matrix](flux-source/hlle-flux-variant-matrix.md).
- NRPy generation/support: `GRHayL/Flux_Source/*.py`,
  `GRHayL/Flux_Source/nrpy/`; route through
  [generated NRPy boundary](flux-source/generated-nrpy-boundary.md).

## Test Paths

- `Unit_Tests/unit_test_HLL_flux.c`
- `Unit_Tests/unit_test_hybrid_flux.c`
- `Unit_Tests/unit_test_tabulated_flux.c`
- `Unit_Tests/unit_test_ET_Legacy_HLL_flux.c`
- `Unit_Tests/unit_test_ET_Legacy_flux_source.c`
- `Unit_Tests/data_gen/unit_test_data_HLL_flux.c`
- `Unit_Tests/data_gen/unit_test_data_hybrid_flux.c`
- `Unit_Tests/data_gen/unit_test_data_tabulated_flux.c`
- `Unit_Tests/data_gen/unit_test_data_ET_Legacy_flux_source.c`

## Key Contracts

- Flux routines expect primitive inputs already reconstructed to faces.
- HLLE variants are split by direction, EOS family, and entropy evolution.
- Source terms require metric derivatives provided by caller-side infrastructure.
- Characteristic speeds are shared by hydrodynamic flux work and induction HLL flux setup.

## Common Edit Routes

- Change characteristic speed: update all affected direction files and tests that cover HLL fluxes and induction flux inputs.
- Add EOS-specific flux variant: add direction implementations, declarations, build-list entries, function-pointer wiring if needed, tests, and docs.
- Change source term: update NRPy source if generated, C output, source tests, and equation map.
- Change reconstruction assumptions: coordinate with `GRHayL/Reconstruction/` and update docs on face orientation.

## Drift Risks

- Direction-specific files can diverge when only one axis is edited.
- Entropy and non-entropy flux variants must stay consistent in conservative field ordering.
- Generated or NRPy-derived expressions can drift from checked-in C if both are not updated together.
- GRHayLib compiles Flux_Source subdirectories directly; new directories require downstream coordination.

## Do Not Duplicate

Keep equations and generated expressions in `docs/raw/derivation.md`, `docs/raw/Flux_Source.dox`, and source. This page is for edit routing.
