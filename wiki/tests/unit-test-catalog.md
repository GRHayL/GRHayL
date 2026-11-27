# Unit Test Catalog

> Compact index of every `Unit_Tests/unit_test_*.c` and `Unit_Tests/compose/*.py` test mapped to module test leaves. · Status: confirmed
> Up: [Tests](index.md)

## Summary

Every standalone C unit-test driver (`Unit_Tests/unit_test_*.c`) and CompOSE Python test (`Unit_Tests/compose/*.py`) is listed below with its subsystem, HDF5/runner mode, and defining module test leaf. Search with `rg -n -i -F '<name>' wiki/tests/unit-test-catalog.md`.

## Detail

### Core, Con2Prim, EOS, Flux, And Induction Tests

| Test | Subsystem & scope | Mode & defining route |
| --- | --- | --- |
| `unit_test_grhayl_core_test_suite.c` | Core EOS bounds, constant atm, 4D metric helpers, ADM det, C2P names. | Default; [Core Tests](../core/tests-and-reference-data.md) |
| `unit_test_compute_conservs_and_Tmunu.c` | Conservative vars and covariant `T_{\mu\nu}` replay. | Default; [Core Tests](../core/tests-and-reference-data.md) |
| `unit_test_code_error.c` | Expected-failure error keys `0`–`90` (`ghl_abort_if_error`). | Default (no-HDF5 skips); [Error Keys](expected-failure-and-error-keys.md) |
| `unit_test_apply_conservative_limits.c` | Con2Prim conservative floors (`tau_fix`, limits). | Default; [Con2Prim Tests](../gems/con2prim/tests-and-reference-data.md) |
| `unit_test_enforce_primitive_limits_and_compute_u0.c` | Primitive floors/ceilings and direct `u0` limiter cases. | Default; [Con2Prim Tests](../gems/con2prim/tests-and-reference-data.md) |
| `unit_test_con2prim_multi_method_hybrid.c` | Hybrid multi-method C2P and `ghl_guess_primitives` replay. | Default; [Con2Prim Tests](../gems/con2prim/tests-and-reference-data.md) |
| `unit_test_hybrid_failure.c` | Hybrid C2P failure codes, backups, `Noble1D_entropy2` checks. | Default, no data; [Con2Prim Tests](../gems/con2prim/tests-and-reference-data.md) |
| `unit_test_con2prim_tabulated.c` | Tabulated C2P replay (NN disabled/enabled, SLy4 + `pyghl`). | HDF5; [Con2Prim Tests](../gems/con2prim/tests-and-reference-data.md) |
| `unit_test_c2p_nn_guess.c` | NN layout, overflow-safe inference, HDF5/fallback checks. | Default + no-HDF5; [Con2Prim Tests](../gems/con2prim/tests-and-reference-data.md) |
| `unit_test_con2prim_debug.c` | Hard-coded tabulated C2P stderr debug dump. | HDF5 manual-only; [Con2Prim Tests](../gems/con2prim/tests-and-reference-data.md) |
| `unit_test_piecewise_polytrope.c` | Piecewise-polytrope `K_ppoly`, `eps_integ_const`, transitions. | Default, no data; [EOS Tests](../gems/eos/tests-and-reference-data.md) |
| `unit_test_tabulated_eos.c` | `simple_table.h5` read, interpolators, bounds, beta-eq. | HDF5; [EOS Tests](../gems/eos/tests-and-reference-data.md) |
| `unit_test_tabulated_eos_compose.c` | Converted CompOSE table EOS/C2P/flux/leakage integration. | HDF5 CI-only; [CompOSE Callback](../gems/neutrinos/compose-eos-leakage-callback.md) |
| `compose/test_compose_to_grhayl.py` | Python CompOSE conversion, regularization, and metadata. | Python CI; [CompOSE Conversion](../gems/neutrinos/compose-eos-conversion-and-regularization.md) |
| `compose/test_compose_failure_coverage.py` | Python CompOSE error-branch tests (CI 100% check). | Python CI; [CompOSE Conversion](../gems/neutrinos/compose-eos-conversion-and-regularization.md) |
| `unit_test_hybrid_flux.c` | Hybrid and hybrid-entropy directional HLLE fluxes. | Default; [Flux_Source Tests](../gems/flux-source/tests-and-reference-data.md) |
| `unit_test_tabulated_flux.c` | Tabulated and tabulated-entropy HLLE fluxes (LS220). | HDF5; [Flux_Source Tests](../gems/flux-source/tests-and-reference-data.md) |
| `unit_test_HLL_flux.c` | Induction vector-potential HLL flux (`B` and `Btilde`). | Default; [Induction Tests](../gems/induction/tests-and-reference-data.md) |
| `unit_test_induction_ccc_ADM.c` | Cell-centered ADM induction interpolation. | Default; [Induction Tests](../gems/induction/tests-and-reference-data.md) |
| `unit_test_induction_ccc_BSSN.c` | Cell-centered BSSN induction interpolation. | Default; [Induction Tests](../gems/induction/tests-and-reference-data.md) |
| `unit_test_induction_vvv_ADM.c` | Vertex-centered ADM induction interpolation. | Default; [Induction Tests](../gems/induction/tests-and-reference-data.md) |

### Neutrinos, Reconstruction, And ET_Legacy Tests

| Test | Subsystem & scope | Mode & defining route |
| --- | --- | --- |
| `unit_test_nrpyleakage_physics.c` | Table-free blocking, Fermi, detailed-balance, stencil checks. | Default + no-HDF5; [Neutrinos Tests](../gems/neutrinos/tests-and-reference-data.md) |
| `unit_test_nrpyleakage_classifier_fallback.c` | Portable finite/Inf/NaN classifier value matrix. | Default + no-HDF5; [Neutrinos Tests](../gems/neutrinos/tests-and-reference-data.md) |
| `unit_test_nrpyleakage_optically_thin_gas.c` | Optically thin leakage evolution + RK lookup fault checks. | HDF5 (SLy4); [Neutrinos Tests](../gems/neutrinos/tests-and-reference-data.md) |
| `unit_test_nrpyleakage_constant_density_sphere.c` | Sphere opacities, optical-depth iteration, neighbor gather. | HDF5 (SLy4); [Neutrinos Tests](../gems/neutrinos/tests-and-reference-data.md) |
| `unit_test_nrpyleakage_luminosities.c` | Fermi branch checks and neutrino luminosity replay. | HDF5 (SLy4); [Neutrinos Tests](../gems/neutrinos/tests-and-reference-data.md) |
| `unit_test_PLM_reconstruction.c` | PLM minmod, MC, and superbee face reconstruction. | Default; [Recon Tests](../gems/reconstruction/tests-and-reference-data.md) |
| `unit_test_WENOZ_reconstruction.c` | WENO-z face reconstruction replay. | Default CI-only; [Recon Tests](../gems/reconstruction/tests-and-reference-data.md) |
| `unit_test_ET_Legacy_conservs.c` | ET_Legacy primitive-to-conservative + `T_{\mu\nu}` replay. | Default; [ET Legacy](et-legacy-comparison.md) |
| `unit_test_ET_Legacy_primitives.c` | ET_Legacy hybrid Con2Prim recovery replay. | Default; [ET Legacy](et-legacy-comparison.md) |
| `unit_test_ET_Legacy_flux_source.c` | ET_Legacy geometric source/flux RHS replay. | Default; [ET Legacy](et-legacy-comparison.md) |
| `unit_test_ET_Legacy_HLL_flux.c` | ET_Legacy Induction vector-potential HLL replay. | Default; [ET Legacy](et-legacy-comparison.md) |
| `unit_test_ET_Legacy_induction_gauge_rhs.c` | ET_Legacy BSSN interpolation and gauge RHS replay. | Default; [ET Legacy](et-legacy-comparison.md) |
| `unit_test_ET_Legacy_reconstruction.c` | ET_Legacy PPM face reconstruction replay. | Default; [ET Legacy](et-legacy-comparison.md) |

Evidence (descriptive behavior): `configure::unit_test_`; corroboration: `Unit_Tests/compose/test_compose_to_grhayl.py` and `Unit_Tests/compose/test_compose_failure_coverage.py`.

## Sources

- [configure](../../configure) - `Unit_Tests/unit_test_*.c` discovery and `--disable-hdf5` filtering
- [.github/run_tests.sh](../../.github/run_tests.sh) - default local test selection
- [Unit_Tests/compose/test_compose_to_grhayl.py](../../Unit_Tests/compose/test_compose_to_grhayl.py) - CompOSE Python conversion tests
- [Unit_Tests/compose/test_compose_failure_coverage.py](../../Unit_Tests/compose/test_compose_failure_coverage.py) - CompOSE Python failure-branch tests

## See Also

- Depends on: [Data Generators And Helpers](data-generators-and-helpers.md)
- Depends on: [Runner And Build Outputs](runner-and-build-outputs.md)
- See also: [Unit Test Coverage And Gap Matrix](unit-test-coverage-and-gap-matrix.md)
