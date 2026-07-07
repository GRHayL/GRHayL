# Con2Prim Gem

## Purpose

Con2Prim recovers primitive variables from conservative variables and provides related primitive-to-conservative helpers. Recovery is numerical and EOS-dependent.

## Read First

- `docs/raw/Con2Prim.dox`
- `docs/raw/derivation.md`
- `wiki/physics/variables-and-conventions.md`
- `wiki/physics/evolution-equation-map.md`

## Public Headers

- `GRHayL/include/ghl_con2prim.h`
- `GRHayL/include/ghl.h`

Key public surface:
- `ghl_con2prim_hybrid_multi_method`
- `ghl_con2prim_tabulated_multi_method`
- `ghl_con2prim_hybrid_select_method`
- `ghl_con2prim_tabulated_select_method`
- `ghl_apply_conservative_limits`
- `ghl_guess_primitives`
- `ghl_enforce_primitive_limits_and_compute_u0`
- `ghl_compute_conservs`
- `ghl_compute_conservs_and_Tmunu`

## Implementation Paths

- Shared routines: `GRHayL/Con2Prim/`
- Hybrid Noble: `GRHayL/Con2Prim/Hybrid/Noble/`
- Hybrid Font1D: `GRHayL/Con2Prim/Hybrid/Font1D/`
- Hybrid Palenzuela1D: `GRHayL/Con2Prim/Hybrid/Palenzuela1D/`
- Tabulated Newman1D: `GRHayL/Con2Prim/Tabulated/Newman1D/`
- Tabulated Noble2D: `GRHayL/Con2Prim/Tabulated/Noble2D/`
- Tabulated Palenzuela1D: `GRHayL/Con2Prim/Tabulated/Palenzuela1D/`
- Tabulated Cerda-Duran path: `GRHayL/Con2Prim/Tabulated/con2prim_CerdaDuran3D.cc`

## Test Paths

- `Unit_Tests/unit_test_con2prim_multi_method_hybrid.c`
- `Unit_Tests/unit_test_con2prim_tabulated.c`
- `Unit_Tests/unit_test_con2prim_debug.c`
- `Unit_Tests/unit_test_hybrid_failure.c`
- `Unit_Tests/unit_test_compute_conservs_and_Tmunu.c`
- `Unit_Tests/data_gen/unit_test_data_con2prim_multi_method_hybrid.c`

## Key Contracts

- Initialize diagnostics with `ghl_initialize_diagnostics`.
- Use `ghl_con2prim_id_t` values in `GRHayL/include/ghl.h` consistently with dispatch and routine-name helpers.
- Respect densitized versus undensitized conservative inputs for each method.
- Always enforce primitive limits and compute `u0` before returning usable primitives.
- Entropy methods must match `params->evolve_entropy` and EOS support.

## Common Edit Routes

- Add hybrid method: update `ghl_con2prim_id_t`, declarations, select/multi-method dispatch, method implementation, build lists, tests, and `docs/raw/Con2Prim.dox`.
- Add tabulated method: same route, plus tabulated EOS variable and bounds assumptions.
- Change conservative limiting: update `apply_conservative_limits`, affected tests, and physics pages if variable meaning changes.
- Change diagnostics: update struct users and tests that assert failure, backup, or iteration behavior.

## Drift Risks

- EOS function-pointer changes can break Con2Prim without local compiler errors in every method.
- Entropy and temperature recovery paths must stay aligned with flux entropy variants and tabulated EOS.
- GRHayLib directly includes Con2Prim public headers and compiles Con2Prim subdirectories; new subdirectories may need downstream build-list coordination.

## Do Not Duplicate

Doxygen pages and source are authority for method references, arguments, and equations. Keep this page as routing and contract guidance only.
