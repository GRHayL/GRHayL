# Workflows

Use these playbooks for cross-module edits. They route to local source of truth and avoid duplicating Doxygen API text.

## EOS Routine

Read first:
- `docs/raw/EOS.dox`
- `docs/raw/GRHayL_Core.dox`
- `GRHayL/include/ghl.h`
- `GRHayL/include/ghl_eos_functions.h`
- `GRHayL/include/ghl_nrpyeos_hybrid.h`
- `GRHayL/include/ghl_nrpyeos_tabulated.h`

Edit paths:
- Hybrid EOS: `GRHayL/EOS/Hybrid/`
- Tabulated EOS: `GRHayL/EOS/Tabulated/`
- Tabulated interpolators: `GRHayL/EOS/Tabulated/interpolators/`
- Table adapters: `GRHayL/EOS/Tabulated/stellarcollapse/`
- Public declarations and function pointers: `GRHayL/include/`
- Build lists: local `make.code.defn` files under changed EOS directories

Tests/data generators:
- `Unit_Tests/unit_test_piecewise_polytrope.c`
- `Unit_Tests/unit_test_tabulated_eos.c`
- `Unit_Tests/test_compute_h_and_cs2.c`
- `Unit_Tests/sample_table/`

Docs to update:
- `docs/raw/EOS.dox` for public behavior or new EOS families
- `docs/raw/GRHayL_Core.dox` if initialization contracts change
- `wiki/gems/eos.md` for routing changes

Pitfalls/contracts:
- Initialize `ghl_eos_parameters` through the core EOS initialization routines so function pointers match the chosen EOS.
- Hybrid and tabulated paths expose different variables and bounds; do not assume a tabulated variable exists in hybrid EOS.
- No-HDF5 builds cannot use HDF5-backed tabulated runtime paths; keep `GHL_DISABLE_HDF5` guards and build lists aligned.
- Changes to pressure, internal energy, entropy, temperature, or `Y_e` can affect Con2Prim and flux paths.

## Con2Prim Method

Read first:
- `wiki/gems/con2prim.md`
- `wiki/gems/con2prim/solver-matrix.md` for solver names, support status, dispatch, and diagnostics
- `wiki/gems/con2prim/recovery-flow.md` for recovery order, backups, and densitized/undensitized boundaries
- `wiki/gems/con2prim/limits-and-conversions.md` for conservative limits, primitive limits, guesses, and conversion helpers
- `wiki/gems/con2prim/tests-and-fixtures.md` for Con2Prim test and fixture routing
- `docs/raw/Con2Prim.dox`
- `docs/raw/derivation.md`
- `GRHayL/include/ghl_con2prim.h`
- `GRHayL/include/ghl.h`

Edit paths:
- Shared pre/post logic: `GRHayL/Con2Prim/`
- Hybrid methods: `GRHayL/Con2Prim/Hybrid/`
- Tabulated methods: `GRHayL/Con2Prim/Tabulated/`
- Public method IDs and diagnostics: `GRHayL/include/ghl.h`, `GRHayL/include/ghl_con2prim.h`
- Build lists: local Con2Prim `make.code.defn` files

Tests/data generators:
- `Unit_Tests/unit_test_con2prim_multi_method_hybrid.c`
- `Unit_Tests/unit_test_con2prim_tabulated.c`
- `Unit_Tests/unit_test_con2prim_debug.c`
- `Unit_Tests/unit_test_hybrid_failure.c`
- `Unit_Tests/data_gen/unit_test_data_con2prim_multi_method_hybrid.c`
- `Unit_Tests/unit_test_compute_conservs_and_Tmunu.c`

Docs to update:
- `docs/raw/Con2Prim.dox` for method availability, EOS support, or references
- `wiki/gems/con2prim.md` for hub routing changes
- `wiki/gems/con2prim/solver-matrix.md` for method support, dispatch, diagnostics, or backup-status changes
- `wiki/gems/con2prim/recovery-flow.md` for driver order, backup behavior, and densitized/undensitized contracts
- `wiki/gems/con2prim/limits-and-conversions.md` for helper contract or conversion changes
- `wiki/gems/con2prim/tests-and-fixtures.md` for fixture, HDF5, or targeted-test changes
- `wiki/physics/variables-and-conventions.md` if variable meaning changes

Pitfalls/contracts:
- `ghl_con2prim_diagnostics` must identify limits, backups, selected routine, and iteration count consistently.
- Some routines take densitized conservatives; others use undensitized values. Check `ghl_undensitize_conservatives` call sites and [recovery flow](gems/con2prim/recovery-flow.md) before editing.
- Entropy methods depend on `params->evolve_entropy` and EOS entropy support.
- Primitive recovery must return bounded primitives and a valid `u0` through `ghl_enforce_primitive_limits_and_compute_u0`.

## Reconstruction Routine

Read first:
- `docs/raw/Reconstruction.dox`
- `GRHayL/include/ghl_reconstruction.h`
- `GRHayL/Reconstruction/`

Edit paths:
- PLM: `GRHayL/Reconstruction/PLM/`
- PPM: `GRHayL/Reconstruction/PPM/`
- WENO-z: `GRHayL/Reconstruction/WENOZ/`
- Public declarations: `GRHayL/include/ghl_reconstruction.h`
- Build lists: local Reconstruction `make.code.defn` files

Tests/data generators:
- `Unit_Tests/unit_test_PLM_reconstruction.c`
- `Unit_Tests/unit_test_WENOZ_reconstruction.c`
- `Unit_Tests/unit_test_ET_Legacy_reconstruction.c`
- `Unit_Tests/data_gen/unit_test_data_PLM_reconstruction.c`
- `Unit_Tests/data_gen/unit_test_data_WENOZ_reconstruction.c`
- `Unit_Tests/data_gen/unit_test_data_ET_Legacy_reconstruction.c`

Docs to update:
- `docs/raw/Reconstruction.dox`
- `wiki/gems/reconstruction.md`
- `wiki/workflows.md` if method routing changes

Pitfalls/contracts:
- Public reconstruction routines return left/right values for the left face of a cell.
- PPM steepening uses pressure, effective Gamma, and parameters in `ghl_parameters`; preserve stencil sizes.
- Reconstruction feeds Flux_Source and Induction face data. Changing output orientation can silently corrupt fluxes.

## Flux/Source Path

Read first:
- `docs/raw/Flux_Source.dox`
- `docs/raw/derivation.md`
- `GRHayL/include/ghl_flux_source.h`
- `GRHayL/Flux_Source/`

Edit paths:
- Characteristic speeds and sources: `GRHayL/Flux_Source/`
- Hybrid HLLE: `GRHayL/Flux_Source/hybrid/`
- Hybrid entropy HLLE: `GRHayL/Flux_Source/hybrid_entropy/`
- Tabulated HLLE: `GRHayL/Flux_Source/tabulated/`
- Tabulated entropy HLLE: `GRHayL/Flux_Source/tabulated_entropy/`
- NRPy generation/support code: `GRHayL/Flux_Source/*.py`, `GRHayL/Flux_Source/nrpy/`
- Public declarations: `GRHayL/include/ghl_flux_source.h`

Tests/data generators:
- `Unit_Tests/unit_test_HLL_flux.c`
- `Unit_Tests/unit_test_hybrid_flux.c`
- `Unit_Tests/unit_test_tabulated_flux.c`
- `Unit_Tests/unit_test_ET_Legacy_HLL_flux.c`
- `Unit_Tests/unit_test_ET_Legacy_flux_source.c`
- `Unit_Tests/data_gen/unit_test_data_HLL_flux.c`
- `Unit_Tests/data_gen/unit_test_data_hybrid_flux.c`
- `Unit_Tests/data_gen/unit_test_data_tabulated_flux.c`

Docs to update:
- `docs/raw/Flux_Source.dox`
- `wiki/gems/flux-source.md`
- `wiki/physics/evolution-equation-map.md` if term mapping changes

Pitfalls/contracts:
- HLLE public entry points are split by direction, EOS family, and entropy evolution.
- Flux inputs are expected at faces; use Reconstruction contracts when changing caller-facing assumptions.
- Source terms require metric derivatives supplied by the surrounding code.
- Characteristic-speed changes also affect Induction HLL flux routing.

## Induction Path

Read first:
- `docs/raw/Induction.dox`
- `docs/raw/derivation.md`
- `GRHayL/include/ghl_induction.h`
- `GRHayL/Induction/`

Edit paths:
- Magnetic HLL fluxes: `GRHayL/Induction/HLL_flux_with_B.c`, `GRHayL/Induction/HLL_flux_with_Btilde.c`
- Gauge/scalar-potential RHS: `GRHayL/Induction/calculate_phitilde_rhs.c`
- Interpolation helpers: `GRHayL/Induction/Interpolators/`
- Public declarations: `GRHayL/include/ghl_induction.h`

Tests/data generators:
- `Unit_Tests/unit_test_induction_ccc_ADM.c`
- `Unit_Tests/unit_test_induction_ccc_BSSN.c`
- `Unit_Tests/unit_test_induction_vvv_ADM.c`
- `Unit_Tests/unit_test_ET_Legacy_induction_gauge_rhs.c`
- `Unit_Tests/compute_A_flux_with_B.c`
- `Unit_Tests/compute_A_flux_with_Btilde.c`
- `Unit_Tests/compute_ccc_ADM.c`
- `Unit_Tests/compute_ccc_BSSN.c`
- `Unit_Tests/compute_vvv_ADM.c`
- `Unit_Tests/data_gen/unit_test_data_induction_interpolation.c`

Docs to update:
- `docs/raw/Induction.dox`
- `wiki/gems/induction.md`
- `wiki/physics/evolution-equation-map.md`

Pitfalls/contracts:
- Magnetic variables use GRHayL's rescaled magnetic-field convention from `docs/raw/derivation.md`.
- `A_i`, `tilde{Phi}`, and staggered magnetic fields have distinct assumed grid locations.
- `ghl_HLL_vars` uses cross-product component ordering, not literal vector index names.
- `ghl_HLL_flux_with_B` densitizes with `psi6`; `ghl_HLL_flux_with_Btilde` expects densitized input.

## NRPyLeakage/Neutrino Path

Read first:
- `wiki/gems/neutrinos.md`
- `wiki/gems/neutrinos/api-and-data.md`
- `wiki/gems/neutrinos/implementation-flow.md`
- `wiki/gems/neutrinos/tests-and-fixtures.md`
- `GRHayL/include/ghl_radiation.h`
- `GRHayL/include/ghl_nrpyleakage.h`
- `GRHayL/Neutrinos/NRPyLeakage/`
- `GRHayL/include/ghl_nrpyeos_tabulated.h`

Edit paths:
- Leakage implementation: `GRHayL/Neutrinos/NRPyLeakage/`
- Radiation structs: `GRHayL/include/ghl_radiation.h`
- Leakage API/constants: `GRHayL/include/ghl_nrpyleakage.h`
- Build lists: `GRHayL/Neutrinos/make.code.defn`, `GRHayL/Neutrinos/NRPyLeakage/make.code.defn`

Tests/data generators:
- `Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c`
- `Unit_Tests/unit_test_nrpyleakage_luminosities.c`
- `Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c`
- `Unit_Tests/nrpyleakage_main.h`

Wiki routes:
- `wiki/gems/neutrinos.md`
- `wiki/gems/neutrinos/api-and-data.md` for radiation structs, public leakage
  entry points, constants ownership, errors, and HDF5/EOS dependency
- `wiki/gems/neutrinos/implementation-flow.md` for Fermi-Dirac,
  PathOfLeastResistance, source terms, luminosities, opacities, and finite
  handling
- `wiki/gems/neutrinos/tests-and-fixtures.md` for `nrpyleakage_main.h`,
  fixture pairs, SLy4 table setup, CI downloads, and key `0`/key `1` behavior

Pitfalls/contracts:
- Leakage routines depend on EOS access to `rho`, `Y_e`, and temperature.
- `ghl_neutrino_opacities` stores electron neutrino, electron antineutrino, and heavy-lepton neutrino entries; optical depths reuse the opacity struct type.
- Units and constants live in `ghl_nrpyleakage.h`; changing them affects all leakage tests.

## Add Unit Test Data

Read first:
- `README.md` CI testing section
- Existing matching `Unit_Tests/unit_test_*.c`
- Existing matching `Unit_Tests/data_gen/unit_test_data_*.c`
- `GRHayL/include/ghl_unit_tests.h`

Edit paths:
- Test executable: `Unit_Tests/unit_test_*.c`
- Data generator: `Unit_Tests/data_gen/unit_test_data_*.c`
- Shared helpers: `Unit_Tests/randomize_*.c`, `Unit_Tests/tabulated_eos_unit_test_helpers.c`
- Sample table inputs: `Unit_Tests/sample_table/`

Tests/data generators:
- Run the changed unit test and its generator when available.
- For perturbation tests, compare output against the generated perturbed bounds, not exact equality unless the existing test does so.

Docs to update:
- Relevant gem page under `wiki/gems/`
- `wiki/test-map.md` is outside this task scope; coordinate with its owner before editing.

Pitfalls/contracts:
- Test data often validates outputs against perturbed input ranges.
- ET_Legacy tests may generate only input data; legacy output comes from downstream validation.
- Do not update shared test maps in parallel-owned files without coordination.

## Doxygen Source Docs

Read first:
- `docs/raw/mainpage.md`
- Relevant `docs/raw/*.dox`
- `Doxyfile`
- Public headers in `GRHayL/include/`

Edit paths:
- Module pages: `docs/raw/*.dox`
- Derivation source: `docs/raw/derivation.md`
- Public API comments: `GRHayL/include/`

Tests/data generators:
- Run Doxygen locally if the build tooling is available.
- At minimum, search links and group names after edits.

Docs to update:
- Keep wiki routing pages compact and point back to Doxygen/source authority.
- Avoid copying long equations or full function docs into wiki pages.

Pitfalls/contracts:
- Doxygen is authoritative for function-level API docs.
- `docs/raw/derivation.md` is authoritative for GRMHD variable definitions and equation concepts.
- Generated output is not source; edit raw docs and headers instead.

## Downstream GRHayLib Impact

Read first:
- `implementations/GRHayLib/README`
- `implementations/GRHayLib/src/GRHayLib.h`
- `implementations/GRHayLib/src/make.code.defn`
- `implementations/GRHayLib/param.ccl`
- `implementations/GRHayLib/schedule.ccl`

Edit paths:
- None in this task scope. Coordinate with the owner before editing `implementations/GRHayLib/`.

Tests/data generators:
- Check whether changed headers are included through `implementations/GRHayLib/src/GRHayLib.h`.
- Check whether new source subdirectories need corresponding downstream `SUBDIRS` entries.

Docs to update:
- Relevant gem page drift-risk section.
- Downstream docs only with coordination.

Pitfalls/contracts:
- GRHayLib directly compiles GRHayL-like subdirectories; new source directories can require downstream build-list changes.
- Public header changes affect Cactus thorns that include `GRHayLib.h`.
- Parameter or EOS initialization changes can require Cactus parameter/schedule updates.
