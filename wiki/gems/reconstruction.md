# Reconstruction Gem

## Purpose

Reconstruction computes shock-aware left/right face values for hydrodynamic variables before flux evaluation.

## Read First

- `docs/raw/Reconstruction.dox`
- `GRHayL/include/ghl_reconstruction.h`
- `wiki/gems/flux-source.md`
- `wiki/gems/induction.md`

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

- PLM: `GRHayL/Reconstruction/PLM/`
- PPM: `GRHayL/Reconstruction/PPM/`
- WENO-z: `GRHayL/Reconstruction/WENOZ/`
- Build lists: `GRHayL/Reconstruction/make.code.defn` and method-local `make.code.defn` files

## Test Paths

- `Unit_Tests/unit_test_PLM_reconstruction.c`
- `Unit_Tests/unit_test_WENOZ_reconstruction.c`
- `Unit_Tests/unit_test_ET_Legacy_reconstruction.c`
- `Unit_Tests/data_gen/unit_test_data_PLM_reconstruction.c`
- `Unit_Tests/data_gen/unit_test_data_WENOZ_reconstruction.c`
- `Unit_Tests/data_gen/unit_test_data_ET_Legacy_reconstruction.c`

## Key Contracts

- Public routines return left/right values for the left face of the current cell.
- PLM methods use compact stencils; PPM and WENO-z use wider stencils.
- PPM steepening depends on pressure, effective Gamma, and PPM parameters in `ghl_parameters`.
- Flux_Source and Induction callers rely on face orientation and stencil size.

## Common Edit Routes

- Add PLM limiter: implement in `GRHayL/Reconstruction/PLM/`, declare it, add tests, and update Doxygen.
- Change PPM behavior: update PPM helpers and both steepened/non-steepened call paths as applicable.
- Change WENO-z behavior: update both main and right-left face variants when the math applies to both.
- Change face orientation: coordinate with Flux_Source, Induction, tests, and Doxygen before landing.

## Drift Risks

- Reconstruction output orientation can break flux tests without obvious local failures.
- PPM parameter changes in `ghl_parameters` affect initialization and downstream GRHayLib parameters.
- New source subdirectories require GRHayLib build-list coordination.

## Do Not Duplicate

Keep method derivation and detailed API docs in `docs/raw/Reconstruction.dox` and source. This page remains edit-routing guidance.
