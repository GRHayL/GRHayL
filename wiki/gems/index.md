# Gem Index

Router for GRHayL core/chalice and gem modules. Use this page to pick docs,
source, headers, tests, and common edit routes before wider search.

## GRHayL_Core

- Purpose: core/chalice layer for shared structs, enums, EOS setup, metric and
  stress-energy helpers, packing/unpacking, error handling, and small utilities.
- Docs path: `docs/raw/GRHayL_Core.dox`, `docs/raw/derivation.md`,
  `docs/raw/mainpage.md`.
- Source path: `GRHayL/GRHayL_Core/`.
- Primary headers: `GRHayL/include/ghl.h`, `GRHayL/include/ghl_io.h`,
  `GRHayL/include/ghl_metric_helpers.h`, `GRHayL/include/ghl_debug.h`,
  `GRHayL/include/ghl_eos_functions.h`,
  `GRHayL/include/ghl_eos_functions_declaration.h`.
- Likely tests: `Unit_Tests/unit_test_grhayl_core_test_suite.c`,
  `Unit_Tests/unit_test_compute_conservs_and_Tmunu.c`, perturbation helpers,
  and randomization helpers.
- Common edit routes: struct or enum changes start in `GRHayL/include/ghl.h`;
  EOS function-pointer dispatch touches core initialization plus EOS/Flux_Source
  implementations; metric helper changes need derivation/core doc review.
- Drift/contract notes: most modules and `implementations/GRHayLib/` consume
  core structs and enums, so public layout and initializer behavior are broad
  compatibility contracts.

## Atmosphere

- Purpose: atmosphere prescriptions for resetting primitive variables to
  atmosphere values.
- KB routes: [hub](atmosphere.md), [prescription contract](atmosphere/prescription-contract.md).
- Docs path: `docs/raw/Atmosphere.dox`.
- Source path: `GRHayL/Atmosphere/`.
- Primary headers: `GRHayL/include/ghl_atmosphere.h`, `GRHayL/include/ghl.h`.
- Likely tests: `Unit_Tests/unit_test_grhayl_core_test_suite.c` directly
  tests constant-atmosphere reset; Con2Prim/EOS limit tests may provide
  indirect coverage.
- Common edit routes: add prescription source, expose it in
  `GRHayL/include/ghl_atmosphere.h`, update `GRHayL/Atmosphere/make.code.defn`,
  then update docs.
- Drift/contract notes: atmosphere values live in `ghl_eos_parameters`; changes
  can affect primitive-limit behavior and downstream setup.

## Con2Prim

- Purpose: conservative-to-primitive inversion methods, method selection,
  diagnostics, primitive/conservative helpers, and limit enforcement.
- KB routes: [hub](con2prim.md), [solver matrix](con2prim/solver-matrix.md),
  [recovery flow](con2prim/recovery-flow.md), [limits and conversions](con2prim/limits-and-conversions.md),
  [tests and fixtures](con2prim/tests-and-fixtures.md).
- Docs path: `docs/raw/Con2Prim.dox`, with physics background in
  `docs/raw/derivation.md`.
- Source path: `GRHayL/Con2Prim/`.
- Primary headers: `GRHayL/include/ghl_con2prim.h`, `GRHayL/include/ghl.h`,
  EOS headers for hybrid/tabulated method support.
- Likely tests: `Unit_Tests/unit_test_con2prim_multi_method_hybrid.c`,
  `Unit_Tests/unit_test_con2prim_tabulated.c`,
  `Unit_Tests/unit_test_apply_conservative_limits.c`,
  `Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c`,
  `Unit_Tests/unit_test_hybrid_failure.c`, related data generators; route
  fixture details through [tests and fixtures](con2prim/tests-and-fixtures.md).
- Common edit routes: new method means add implementation, expose/select it in
  `GRHayL/include/ghl_con2prim.h` and `ghl_con2prim_id_t`, update method-name
  mapping, update docs/tests, and check `GRHayLib` parameter mapping. Use
  [solver matrix](con2prim/solver-matrix.md) for method-support evidence.
- Drift/contract notes: solver IDs, backup routine behavior, entropy/temp
  assumptions, and tabulated EOS support are externally visible through
  downstream configuration. Use [recovery flow](con2prim/recovery-flow.md) for
  diagnostics, dispatch, backups, and densitized/undensitized boundaries.

## EOS

- Purpose: simple, hybrid, and tabulated equation-of-state support, including
  bounds, table reading, function-pointer dispatch, enthalpy, and sound speed.
- KB routes: [hub](eos.md), [initialization and dispatch](eos/initialization-and-dispatch.md),
  [hybrid piecewise-polytrope](eos/hybrid-piecewise-polytrope.md),
  [tabulated table contract](eos/tabulated-table-contract.md),
  [tabulated interpolation and bounds](eos/tabulated-interpolation-and-bounds.md),
  [tests and fixtures](eos/tests-and-fixtures.md).
- Docs path: `docs/raw/EOS.dox`, plus HDF5 build notes in `README.md`.
- Source path: `GRHayL/EOS/`.
- Primary headers: `GRHayL/include/ghl.h`,
  `GRHayL/include/ghl_eos_functions.h`,
  `GRHayL/include/ghl_eos_functions_declaration.h`,
  `GRHayL/include/ghl_nrpyeos_hybrid.h`,
  `GRHayL/include/ghl_nrpyeos_tabulated.h`.
- Likely tests: `Unit_Tests/unit_test_piecewise_polytrope.c`,
  `Unit_Tests/unit_test_tabulated_eos.c`,
  `Unit_Tests/test_compute_h_and_cs2.c`,
  `Unit_Tests/tabulated_eos_unit_test_helpers.c`, `Unit_Tests/sample_table/`.
- Common edit routes: hybrid changes touch `GRHayL/EOS/Hybrid/`, hybrid header,
  and piecewise-polytrope tests; tabulated changes touch `GRHayL/EOS/Tabulated/`,
  HDF5 build setup, tabulated tests, and Neutrinos/Con2Prim callers.
- Drift/contract notes: `GHL_DISABLE_HDF5`, omitted tabulated sources, and
  function-pointer initialization must stay consistent across build scripts,
  docs, and downstream `GRHayLib` HDF5 requirements.

## Flux_Source

- Purpose: characteristic speeds, GRMHD source terms, and EOS-specific HLLE
  fluxes for hydrodynamic evolution.
- Docs path: `docs/raw/Flux_Source.dox`, with equations in
  `docs/raw/derivation.md`.
- Source path: `GRHayL/Flux_Source/`.
- Primary headers: `GRHayL/include/ghl_flux_source.h`,
  `GRHayL/include/ghl_eos_functions.h`, `GRHayL/include/ghl.h`.
- Likely tests: `Unit_Tests/unit_test_hybrid_flux.c`,
  `Unit_Tests/unit_test_tabulated_flux.c`,
  `Unit_Tests/unit_test_ET_Legacy_flux_source.c`,
  `Unit_Tests/unit_test_ET_Legacy_HLL_flux.c`, flux data generators.
- Common edit routes: change generated equations or C kernels together; update
  EOS-specific flux variants, function-pointer setup, tests, and Doxygen.
- Drift/contract notes: flux inputs expect reconstructed face primitives.
  Entropy and tabulated variants must match EOS initialization choices.

## Induction

- Purpose: vector-potential magnetic fluxes, scalar-potential RHS, Lorenz gauge
  terms, and interpolation between metric/vector-potential staggerings.
- Docs path: `docs/raw/Induction.dox`, with magnetic-field conventions in
  `docs/raw/derivation.md`.
- Source path: `GRHayL/Induction/`.
- Primary headers: `GRHayL/include/ghl_induction.h`,
  `GRHayL/include/ghl_flux_source.h`, `GRHayL/include/ghl.h`.
- Likely tests: `Unit_Tests/unit_test_HLL_flux.c`,
  `Unit_Tests/unit_test_induction_ccc_ADM.c`,
  `Unit_Tests/unit_test_induction_ccc_BSSN.c`,
  `Unit_Tests/unit_test_induction_vvv_ADM.c`,
  `Unit_Tests/unit_test_ET_Legacy_induction_gauge_rhs.c`.
- Common edit routes: interpolation changes start in
  `GRHayL/Induction/Interpolators/` and `GRHayL/include/ghl_induction.h`;
  magnetic flux changes need characteristic-speed and HLL tests.
- Drift/contract notes: staggered-grid assumptions are public caller contracts.
  Keep docs, stencil shapes, and tests synchronized.

## Neutrinos

- Purpose: NRPyLeakage neutrino opacities, luminosities, optical-depth updates,
  and source terms.
- KB routes: [hub](neutrinos.md), [API and data](neutrinos/api-and-data.md),
  [implementation flow](neutrinos/implementation-flow.md), and
  [tests and fixtures](neutrinos/tests-and-fixtures.md).
- Docs path: coverage gap: no dedicated `docs/raw/Neutrinos.dox` exists;
  `GRHayL/include/ghl_radiation.h` defines the Doxygen group.
- Source path: `GRHayL/Neutrinos/`.
- Primary headers: `GRHayL/include/ghl_radiation.h`,
  `GRHayL/include/ghl_nrpyleakage.h`,
  `GRHayL/include/ghl_nrpyeos_tabulated.h`.
- Likely tests: `Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c`,
  `Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c`,
  `Unit_Tests/unit_test_nrpyleakage_luminosities.c`,
  `Unit_Tests/nrpyleakage_main.h`; route fixture details through
  [tests and fixtures](neutrinos/tests-and-fixtures.md).
- Common edit routes: add or change leakage routines in
  `GRHayL/Neutrinos/NRPyLeakage/`, expose API in leakage/radiation headers,
  update tests and downstream header aggregation. Use
  [API and data](neutrinos/api-and-data.md) for public structs, entry points,
  and HDF5/EOS contracts, and [implementation flow](neutrinos/implementation-flow.md)
  for the five built source files.
- Drift/contract notes: leakage uses tabulated EOS quantities and table-backed
  test data. HDF5/EOS changes can break Neutrinos even if leakage source is
  untouched.

## Reconstruction

- Purpose: PLM, PPM, and WENO-z shock-capturing reconstruction routines that
  produce face values for flux calculations.
- Docs path: `docs/raw/Reconstruction.dox`.
- Source path: `GRHayL/Reconstruction/`.
- Primary headers: `GRHayL/include/ghl_reconstruction.h`,
  `GRHayL/include/ghl.h`.
- Likely tests: `Unit_Tests/unit_test_PLM_reconstruction.c`,
  `Unit_Tests/unit_test_WENOZ_reconstruction.c`,
  `Unit_Tests/unit_test_ET_Legacy_reconstruction.c`, reconstruction data
  generators.
- Common edit routes: PLM changes go through `GRHayL/Reconstruction/PLM/`,
  PPM through `GRHayL/Reconstruction/PPM/`, WENO-z through
  `GRHayL/Reconstruction/WENOZ/`; update header, docs, tests, and
  `make.code.defn`.
- Drift/contract notes: PPM parameters live in `ghl_parameters`; stencil sizes,
  face orientation, and left/right naming are external caller contracts.
  Coverage gap: no dedicated PPM unit test file is obvious in `Unit_Tests/`.
