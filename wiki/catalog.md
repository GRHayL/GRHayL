# GRHayL Query Catalog

Use this catalog to route common terms and aliases to KB pages, source or
headers, and docs or tests. It is a router, not a replacement for Doxygen or
source.

| Term or alias | Route to KB page | Source or header path | Doc or test path |
| --- | --- | --- | --- |
| core, chalice, structs, packing, unpacking | [Public API Map](public-api-map.md) | `GRHayL/GRHayL_Core/`, `GRHayL/include/ghl.h`, `GRHayL/include/ghl_metric_helpers.h` | `docs/raw/GRHayL_Core.dox`, `Unit_Tests/unit_test_grhayl_core_test_suite.c` |
| EOS initialization, function pointers | [Public API Map](public-api-map.md) | `GRHayL/GRHayL_Core/initialize_eos.c`, `GRHayL/include/ghl.h`, `GRHayL/include/ghl_eos_functions.h` | `docs/raw/GRHayL_Core.dox`, `docs/raw/EOS.dox` |
| simple EOS, ideal gas, ideal fluid | [Gems](gems/index.md) | `GRHayL/EOS/`, `GRHayL/GRHayL_Core/initialize_eos.c`, `GRHayL/include/ghl.h` | `docs/raw/EOS.dox`, `Unit_Tests/test_compute_h_and_cs2.c` |
| hybrid EOS, piecewise polytrope, cold pressure | [Gems](gems/index.md) | `GRHayL/EOS/Hybrid/`, `GRHayL/include/ghl_nrpyeos_hybrid.h`, `GRHayL/include/ghl_eos_functions.h` | `docs/raw/EOS.dox`, `Unit_Tests/unit_test_piecewise_polytrope.c` |
| tabulated EOS, HDF5 table, table bounds | [Gems](gems/index.md) | `GRHayL/EOS/Tabulated/`, `GRHayL/include/ghl_nrpyeos_tabulated.h`, `GRHayL/include/ghl_eos_functions.h` | `docs/raw/EOS.dox`, `Unit_Tests/unit_test_tabulated_eos.c`, `Unit_Tests/sample_table/` |
| HDF5 off, no HDF5, `GHL_DISABLE_HDF5`, `--disable-hdf5` | [Build And CI](build-and-ci.md) | `configure`, `GRHayL/EOS/Tabulated/`, `GRHayL/Con2Prim/Tabulated/`, `GRHayL/Flux_Source/tabulated/` | `README.md`, `.github/run_tests.sh` |
| c2p, primitive recovery, conservative-to-primitive | [Gems](gems/index.md) | `GRHayL/Con2Prim/`, `GRHayL/include/ghl_con2prim.h` | `docs/raw/Con2Prim.dox`, `Unit_Tests/unit_test_con2prim_multi_method_hybrid.c`, `Unit_Tests/unit_test_con2prim_tabulated.c` |
| Noble, Font, Palenzuela, Newman solvers | [Gems](gems/index.md) | `GRHayL/Con2Prim/Hybrid/`, `GRHayL/Con2Prim/Tabulated/`, `GRHayL/Con2Prim/utils_Noble.h`, `GRHayL/Con2Prim/utils_Palenzuela1D.h` | `docs/raw/Con2Prim.dox`, `Unit_Tests/unit_test_con2prim_multi_method_hybrid.c`, `Unit_Tests/unit_test_con2prim_tabulated.c` |
| conservative limits, primitive limits, speed limiter | [Gems](gems/index.md) | `GRHayL/Con2Prim/apply_conservative_limits.c`, `GRHayL/Con2Prim/enforce_primitive_limits_and_compute_u0.c`, `GRHayL/include/ghl_con2prim.h` | `Unit_Tests/unit_test_apply_conservative_limits.c`, `Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c` |
| atmosphere, constant atmosphere | [Gems](gems/index.md) | `GRHayL/Atmosphere/`, `GRHayL/include/ghl_atmosphere.h` | `docs/raw/Atmosphere.dox`, `Unit_Tests/unit_test_grhayl_core_test_suite.c` |
| HLLE flux, HLL flux, characteristic speeds | [Evolution Equation Map](physics/evolution-equation-map.md) | `GRHayL/Flux_Source/`, `GRHayL/Flux_Source/hybrid/`, `GRHayL/Flux_Source/tabulated/`, `GRHayL/include/ghl_flux_source.h` | `docs/raw/Flux_Source.dox`, `Unit_Tests/unit_test_hybrid_flux.c`, `Unit_Tests/unit_test_tabulated_flux.c`, `Unit_Tests/unit_test_HLL_flux.c` |
| source terms, RHS, metric derivatives | [Evolution Equation Map](physics/evolution-equation-map.md) | `GRHayL/Flux_Source/ghl_calculate_source_terms.c`, `GRHayL/Flux_Source/GRHayL_rhs.py`, `GRHayL/include/ghl_flux_source.h` | `docs/raw/Flux_Source.dox`, `docs/raw/derivation.md`, `Unit_Tests/unit_test_ET_Legacy_flux_source.c` |
| generated flux/source code, NRPy, Python generators | [Generated Boundaries](generated-boundaries.md) | `GRHayL/Flux_Source/GRHayL_rhs.py`, `GRHayL/Flux_Source/GRMHD_equations_new_version.py`, `GRHayL/Flux_Source/nrpy/` | `Doxyfile`, `docs/raw/Flux_Source.dox` |
| induction, vector potential, `A_i`, `phitilde`, `tildePhi` | [Evolution Equation Map](physics/evolution-equation-map.md) | `GRHayL/Induction/`, `GRHayL/include/ghl_induction.h` | `docs/raw/Induction.dox`, `Unit_Tests/unit_test_induction_ccc_ADM.c`, `Unit_Tests/unit_test_induction_vvv_ADM.c`, `Unit_Tests/unit_test_induction_ccc_BSSN.c` |
| induction interpolation, ADM, BSSN, staggered grids | [Evolution Equation Map](physics/evolution-equation-map.md) | `GRHayL/Induction/Interpolators/`, `GRHayL/include/ghl_induction.h` | `docs/raw/Induction.dox`, `Unit_Tests/unit_test_induction_ccc_ADM.c`, `Unit_Tests/unit_test_induction_vvv_ADM.c`, `Unit_Tests/unit_test_induction_ccc_BSSN.c` |
| NRPyLeakage, neutrinos, opacities, luminosities, optical depths | [Gems](gems/index.md) | `GRHayL/Neutrinos/NRPyLeakage/`, `GRHayL/include/ghl_radiation.h`, `GRHayL/include/ghl_nrpyleakage.h` | `Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c`, `Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c`, `Unit_Tests/unit_test_nrpyleakage_luminosities.c` |
| WENOZ, WENO-z, reconstruction | [Gems](gems/index.md) | `GRHayL/Reconstruction/WENOZ/`, `GRHayL/include/ghl_reconstruction.h` | `docs/raw/Reconstruction.dox`, `Unit_Tests/unit_test_WENOZ_reconstruction.c` |
| PPM, piecewise parabolic method, steepening, flattening | [Gems](gems/index.md) | `GRHayL/Reconstruction/PPM/`, `GRHayL/include/ghl_reconstruction.h` | `docs/raw/Reconstruction.dox` |
| PLM, minmod, maxmod, monotonized central, superbee | [Gems](gems/index.md) | `GRHayL/Reconstruction/PLM/`, `GRHayL/include/ghl_reconstruction.h` | `docs/raw/Reconstruction.dox`, `Unit_Tests/unit_test_PLM_reconstruction.c` |
| primitive variables, conservative variables, magnetic rescaling | [Physics Variables](physics/variables-and-conventions.md) | `GRHayL/include/ghl.h` | `docs/raw/derivation.md`, `docs/raw/mainpage.md` |
| stress-energy tensor, `Tmunu`, small b, ADM auxiliaries | [Public API Map](public-api-map.md) | `GRHayL/GRHayL_Core/compute_TDNmunu.c`, `GRHayL/GRHayL_Core/compute_TUPmunu.c`, `GRHayL/GRHayL_Core/compute_smallb_and_b2.c`, `GRHayL/include/ghl.h` | `docs/raw/GRHayL_Core.dox`, `docs/raw/derivation.md`, `Unit_Tests/unit_test_compute_conservs_and_Tmunu.c` |
| metric helpers, raise, lower, vector square | [Public API Map](public-api-map.md) | `GRHayL/include/ghl_metric_helpers.h`, `GRHayL/GRHayL_Core/compute_ADM_auxiliaries.c` | `docs/raw/GRHayL_Core.dox`, `Unit_Tests/randomize_metric.c` |
| ET_Legacy, IllinoisGRMHD, Einstein Toolkit | [Test Map](test-map.md) | `Unit_Tests/unit_test_ET_Legacy_*.c`, `Unit_Tests/data_gen/unit_test_data_ET_Legacy_*.c`, `implementations/GRHayLib/` | `README.md`, `docs/raw/mainpage.md`, `.github/run_tests.sh` |
| GRHayLib, downstream implementation, Cactus thorn files | [Source Map](source-map.md) | `implementations/GRHayLib/`, `implementations/GRHayLib/src/` | `implementations/GRHayLib/README`, `implementations/GRHayLib/doc/documentation.tex` |
| tests, unit tests, perturbation checks | [Test Map](test-map.md) | `Unit_Tests/`, `GRHayL/include/ghl_unit_tests.h` | `README.md`, `docs/raw/mainpage.md`, `.github/run_tests.sh` |
| fixtures, reference data, TestData/Test_Data drift | [Contradictions](contradictions.md) | `Unit_Tests/data_gen/`, `.github/run_tests.sh` | `README.md`, `docs/raw/mainpage.md`, `Unit_Tests/sample_table/` |
| sample EOS table, generated table | [Test Map](test-map.md) | `Unit_Tests/sample_table/generate_simple_table.py` | `Unit_Tests/sample_table/table_info.txt`, `Unit_Tests/sample_table/` |
| Doxygen, docs output, generated docs boundary | [Generated Boundaries](generated-boundaries.md) | `Doxyfile`, `docs/raw/` | `docs/raw/mainpage.md`, `docs/raw/*.dox` |
| `make.code.defn`, build source lists, generated Makefile inputs | [Generated Boundaries](generated-boundaries.md) | `GRHayL/make.code.defn`, `GRHayL/*/make.code.defn`, `GRHayL/include/make.code.defn` | `configure`, `generate_makefile.sh` |
| CI workflows, GitHub Actions, test runner | [Build And CI](build-and-ci.md) | `.github/workflows/`, `.github/actions/`, `.github/run_tests.sh` | `README.md`, `Unit_Tests/` |
| licensing, citations, inherited sources | [Source Map](source-map.md) | `LICENSE`, `AUTHORS` | `docs/raw/license.md`, `docs/raw/ref.bib`, `docs/raw/mainpage.md` |

## Catalog Rules

- Add aliases here when a reader may not know the owning gem.
- Route to KB pages first, then source/header, then docs/tests.
- Keep paths repo-relative.
- Do not copy Doxygen tables or source bodies; link to them.
