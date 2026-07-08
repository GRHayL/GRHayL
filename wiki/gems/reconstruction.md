# Reconstruction Gem

## Purpose

Reconstruction routes shock-aware left/right face-value methods used before
flux evaluation. Keep method detail in the child pages, source, and read-only
Doxygen evidence.

## Read First

- [Face and stencil contract](reconstruction/face-and-stencil-contract.md)
- [PLM limiters](reconstruction/plm-limiters.md) for minmod, maxmod, MC, and
  superbee lookup.
- [PPM flow](reconstruction/ppm-flow.md)
- [WENOZ contract](reconstruction/wenoz-contract.md)
- [Tests and fixtures](reconstruction/tests-and-fixtures.md)
- `docs/raw/Reconstruction.dox` is read-only evidence for this KB work.
- `GRHayL/include/ghl_reconstruction.h`
- `GRHayL/Reconstruction/`

## Public Headers

- `GRHayL/include/ghl_reconstruction.h`
- `GRHayL/include/ghl.h`

Key public surface:
- `ghl_minmod_reconstruction`
- `ghl_mc_reconstruction`
- `ghl_superbee_reconstruction`
- `ghl_ppm_reconstruction`
- `ghl_ppm_reconstruction_with_steepening`
- `ghl_wenoz_reconstruction`
- `ghl_wenoz_reconstruction_right_left_faces`

## Implementation Paths

- PLM: `GRHayL/Reconstruction/PLM/`; route through [PLM limiters](reconstruction/plm-limiters.md).
- PPM: `GRHayL/Reconstruction/PPM/`; route through [PPM flow](reconstruction/ppm-flow.md).
- WENO-z: `GRHayL/Reconstruction/WENOZ/`; route through [WENOZ contract](reconstruction/wenoz-contract.md).
- Build lists: `GRHayL/Reconstruction/make.code.defn` and method-local `make.code.defn` files

## Test Paths

- `Unit_Tests/unit_test_PLM_reconstruction.c`
- `Unit_Tests/unit_test_WENOZ_reconstruction.c`
- `Unit_Tests/unit_test_ET_Legacy_reconstruction.c`
- `Unit_Tests/data_gen/unit_test_data_PLM_reconstruction.c`
- `Unit_Tests/data_gen/unit_test_data_WENOZ_reconstruction.c`
- `Unit_Tests/data_gen/unit_test_data_ET_Legacy_reconstruction.c`
- Route fixture and CI details through [Tests and fixtures](reconstruction/tests-and-fixtures.md).

## Key Contracts

- Public routines return left/right values for the left face of the current cell; see [Face and stencil contract](reconstruction/face-and-stencil-contract.md).
- PLM, PPM, and WENO-z have different stencil widths; route by method page before changing callers.
- PPM steepening depends on pressure, effective Gamma, and PPM parameters in `ghl_parameters`; see [PPM flow](reconstruction/ppm-flow.md).
- Flux_Source and Induction callers rely on face orientation and stencil size.

## Common Edit Routes

- Add PLM limiter: start with [PLM limiters](reconstruction/plm-limiters.md), then implement in `GRHayL/Reconstruction/PLM/`, declare it, and update tests/routes as needed.
- Change PPM behavior: update PPM helpers and both steepened/non-steepened call paths as applicable.
- Change WENO-z behavior: update both main and right-left face variants when the math applies to both.
- Change face orientation: start with [Face and stencil contract](reconstruction/face-and-stencil-contract.md), then coordinate Flux_Source, Induction, and tests before landing.

## Drift Risks

- Reconstruction output orientation can break flux tests without obvious local failures.
- PPM parameter changes in `ghl_parameters` affect initialization and downstream GRHayLib parameters.
- New source subdirectories require GRHayLib build-list coordination.

## Do Not Duplicate

Keep method derivation and detailed API reference in `docs/raw/Reconstruction.dox` and source. This page remains routing guidance.
