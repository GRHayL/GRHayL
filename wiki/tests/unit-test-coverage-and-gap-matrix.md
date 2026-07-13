# Unit Test Coverage And Gap Matrix

This page maps Unit_Tests evidence to behavior families and known weak areas.
It is not Doxygen, not API documentation, and not a replacement for source,
runner scripts, or per-gem test pages.

Use [Test Map](../test-map.md) for the full file inventory. Use the
cross-cutting Unit_Tests pages for fixture/run details:
[runner and generated artifacts](runner-and-generated-artifacts.md),
[fixture lifecycle and harness contract](fixture-lifecycle-and-harness-contract.md),
[HDF5 sample tables](hdf5-sample-tables.md),
[ET Legacy comparison contract](et-legacy-comparison-contract.md), and
[expected-failure and error keys](expected-failure-and-error-keys.md).

Test helpers are test-only evidence. Helper APIs and helper files under
`Unit_Tests/` are not production public API. `ghl_unit_tests.h` is nevertheless
listed by the install manifest; this creates installed-header presence, not a
linkable production helper surface.

Matrix evidence types describe repository structure: `direct test` means the
test body directly calls/exercises the named family when executed;
`workflow-only` means YAML selects a command. Neither label records a passing
historical run. Missing direct evidence is a `coverage-gap`, not proof of a
behavior defect.

## Coverage Matrix

| Behavior family | Unit_Tests evidence | Evidence type | Main route | Conservative reading |
| --- | --- | --- | --- | --- |
| Core/chalice helpers | `unit_test_grhayl_core_test_suite.c`; indirect use in `unit_test_compute_conservs_and_Tmunu.c`, `unit_test_ET_Legacy_conservs.c`, Con2Prim and flux tests | direct test; replay fixture test; generator evidence | [Core tests and fixtures](../core/tests-and-fixtures.md) | Direct Core coverage includes EOS setup defaults, 4D metric helpers, ADM metric init, atmosphere reset for simple/hybrid EOS, and Con2Prim routine names. It does not prove every Core helper. |
| EOS | `unit_test_piecewise_polytrope.c`, `unit_test_tabulated_eos.c`, `unit_test_code_error.c`, `test_compute_h_and_cs2.c`, `tabulated_eos_unit_test_helpers.c` | direct test; helper-only; expected-failure; HDF5-only; no-HDF5 skip | [EOS tests and fixtures](../gems/eos/tests-and-fixtures.md), [HDF5 sample tables](hdf5-sample-tables.md) | Hybrid piecewise-polytrope and tabulated interpolation/table paths have direct tests. Tabulated atmosphere default reset is still weak because the Core test has a commented/TODO branch. |
| Con2Prim | `unit_test_con2prim_multi_method_hybrid.c`, `unit_test_con2prim_tabulated.c`, `unit_test_c2p_nn_guess.c`, `unit_test_hybrid_failure.c`, `unit_test_apply_conservative_limits.c`, `unit_test_enforce_primitive_limits_and_compute_u0.c`, `unit_test_compute_conservs_and_Tmunu.c`, `unit_test_con2prim_debug.c`, ET Legacy conservs/primitives | direct test; replay fixture test; expected-failure; HDF5-only; manual-only | [Con2Prim tests and fixtures](../gems/con2prim/tests-and-fixtures.md) | Hybrid, tabulated, NN primitive-guess, limit, failure, and stress-energy/conservative paths are covered through selected direct tests and fixtures. `unit_test_con2prim_debug.c` is manual-only, not default runner coverage. |
| Flux_Source | `unit_test_hybrid_flux.c`, `unit_test_tabulated_flux.c`, `unit_test_ET_Legacy_flux_source.c`, `data_gen/unit_test_data_hybrid_flux.c`, `data_gen/unit_test_data_tabulated_flux.c`, `data_gen/unit_test_data_ET_Legacy_flux_source.c` | direct test; replay fixture test; generator evidence; HDF5-only | [Flux Source hub](../gems/flux-source.md), [Flux_Source tests and fixtures](../gems/flux-source/tests-and-fixtures.md), [HLLE flux variant matrix](../gems/flux-source/hlle-flux-variant-matrix.md) | Hybrid/tabulated direct variant HLLE functions and ET Legacy flux/source replay are covered. Generic HLLE dispatch/function-pointer surface is not proven by these Unit_Tests direct variant tests. `unit_test_HLL_flux.c` is Induction vector-potential HLL evidence, not Flux_Source HLLE coverage. |
| Induction | `unit_test_HLL_flux.c`, `unit_test_ET_Legacy_HLL_flux.c`, `unit_test_induction_ccc_ADM.c`, `unit_test_induction_ccc_BSSN.c`, `unit_test_induction_vvv_ADM.c`, `unit_test_ET_Legacy_induction_gauge_rhs.c`, top-level `compute_*` helpers | direct test; replay fixture test; helper-only; generator evidence | [Induction tests and fixtures](../gems/induction/tests-and-fixtures.md), [Induction verification workflows](../gems/induction/verification-workflows.md) | Vector-potential HLL, selected interpolation, and ET Legacy gauge RHS paths are covered. No visible standalone `phitilde_rhs`-only test and no visible `vvv_BSSN` Induction test appear in current evidence. |
| Neutrinos | `unit_test_nrpyleakage_optically_thin_gas.c`, `unit_test_nrpyleakage_constant_density_sphere.c`, `unit_test_nrpyleakage_luminosities.c`, `nrpyleakage_main.h`, `unit_test_code_error.c` | replay fixture test; helper-only; expected-failure; HDF5-only; no-HDF5 skip | [Neutrinos tests and fixtures](../gems/neutrinos/tests-and-fixtures.md), [HDF5 sample tables](hdf5-sample-tables.md) | NRPyLeakage replay uses shared key `1`; local key `0` generation is not normal CI replay. Invalid Fermi-Dirac key paths route through expected-failure coverage. |
| Reconstruction | `unit_test_PLM_reconstruction.c`, `unit_test_WENOZ_reconstruction.c`, `unit_test_ET_Legacy_reconstruction.c`, matching data generators | direct test; replay fixture test; generator evidence; workflow-only | [Reconstruction tests and fixtures](../gems/reconstruction/tests-and-fixtures.md) | PLM and WENOZ have direct fixture tests, but `run_tests.sh` runs PLM and ET Legacy while workflows include WENOZ. ET Legacy exercises PPM paths; no dedicated `unit_test_PPM_reconstruction.c` or PPM data generator is visible. |
| ET Legacy | `unit_test_ET_Legacy_conservs.c`, `unit_test_ET_Legacy_primitives.c`, `unit_test_ET_Legacy_flux_source.c`, `unit_test_ET_Legacy_HLL_flux.c`, `unit_test_ET_Legacy_induction_gauge_rhs.c`, `unit_test_ET_Legacy_reconstruction.c`, matching data generators | replay fixture test; generator evidence | [ET Legacy comparison contract](et-legacy-comparison-contract.md), [Test Map](../test-map.md) | ET Legacy tests are upstream GRHayL regression evidence. Visible local generators often write input-side data; normal runs download trusted outputs. Do not treat them as direct GRHayLib/Cactus verification. |
| Error paths | `unit_test_code_error.c` keys `0..85`; runner loop in `.github/run_tests.sh`; workflows' expected-failure jobs | expected-failure; HDF5-only; no-HDF5 skip; workflow-only | [expected-failure and error keys](expected-failure-and-error-keys.md), [Core tests and fixtures](../core/tests-and-fixtures.md), [EOS tests and fixtures](../gems/eos/tests-and-fixtures.md) | These tests intentionally expect failure/error returns, including NN loader/init keys `83..85`. In no-HDF5 builds, HDF5-only keys are skipped as expected by the harness, not counted as ordinary behavior coverage. |
| Fixtures/helpers | `ghl_unit_tests.h`; `pert_test_fail_conservatives.c`, `pert_test_fail_primitives.c`, `pert_test_fail_stress_energy.c`, `randomize_metric.c`, `randomize_primitives.c`, `compute_A_flux_with_B*.c`, `compute_ccc_*.c`, `compute_vvv_ADM.c`, EOS helper files | installed header; helper-only; generator evidence | [fixture lifecycle and harness contract](fixture-lifecycle-and-harness-contract.md), [Build And CI](../build-and-ci.md), per-gem test pages | Header's inline helpers are caller-compiled. Non-inline helper definitions live outside library manifests, and some declared binary read/write helpers plus `ghl_initial_random_data` have no visible definition. Installed presence must not be cited as linkable production API. |
| Data generators | `Unit_Tests/data_gen/unit_test_data_*.c`; generation modes in `unit_test_con2prim_tabulated.c` and NRPyLeakage harness | generator evidence; HDF5-only where table-backed | [runner and generated artifacts](runner-and-generated-artifacts.md), [fixture lifecycle and harness contract](fixture-lifecycle-and-harness-contract.md) | Generator presence does not mean CI regenerates trusted data. `run_tests.sh` downloads many root-level `.bin` fixtures before replay. |
| Sample tables | `Unit_Tests/sample_table/Hempel_*.h5`, `generate_simple_table.py`, `table_info.txt`; runner downloads `simple_table.h5`, LS220, and SLy4 tables | sample-table asset; HDF5-only | [HDF5 sample tables](hdf5-sample-tables.md), [EOS tests and fixtures](../gems/eos/tests-and-fixtures.md) | Checked-in sample-table assets and downloaded runtime tables support EOS, Flux_Source, Con2Prim, Neutrinos, and error-path tests. Table API behavior belongs in EOS pages. |
| Run modes | `configure`, `make tests`, `make datagen`, `.github/run_tests.sh`, `.github/workflows/*.yml` | workflow-only; manual-only; HDF5-only; no-HDF5 skip | [runner and generated artifacts](runner-and-generated-artifacts.md), [Test Map](../test-map.md) | `configure --disable-hdf5` filters tabulated/HDF5 targets. `run_tests.sh` is broad replay coverage but omits some workflow matrix tests such as WENOZ. Manual-only debug routes must not be counted as default runner coverage. |

## Weak And Gap Callouts

- Incomplete direct Core pack/unpack round-trip coverage: visible Unit_Tests do
  not directly round-trip all `ghl_initialize_*`/`ghl_return_*` pairs. Current
  evidence includes indirect repacking in ET Legacy conservs and broader Core
  initialization use.
- Flux_Source generic HLLE dispatch/function-pointer surface is not proven by
  Unit_Tests direct variant tests. Existing tests select concrete hybrid,
  hybrid-entropy, tabulated, or tabulated-entropy direction functions.
- Helper APIs/helpers in `Unit_Tests/` are test-only and not production public
  API. Do not cite helper-only files as library API coverage. Installed
  `ghl_unit_tests.h` non-inline declarations are a `coverage-gap`/surface
  contradiction, not proof of implementation behavior.
- Core weak areas already routed in [Core tests and fixtures](../core/tests-and-fixtures.md):
  clamp helpers, standalone 3D metric helpers, direct `ghl_compute_TUPmunu`
  fixture validation, and tabulated atmosphere reset branch.
- Reconstruction weak area already routed in
  [Reconstruction tests and fixtures](../gems/reconstruction/tests-and-fixtures.md):
  no dedicated PPM unit-test file or PPM data generator is visible, although ET
  Legacy reconstruction exercises PPM paths.
- Induction weak areas already routed in
  [Induction tests and fixtures](../gems/induction/tests-and-fixtures.md): no
  visible standalone `phitilde_rhs`-only test and no visible `vvv_BSSN`
  interpolation test.

## Repo-Local References

- [wiki/test-map.md](../test-map.md)
- [wiki/tests/index.md](index.md)
- [wiki/tests/runner-and-generated-artifacts.md](runner-and-generated-artifacts.md)
- [wiki/tests/fixture-lifecycle-and-harness-contract.md](fixture-lifecycle-and-harness-contract.md)
- [wiki/tests/hdf5-sample-tables.md](hdf5-sample-tables.md)
- [wiki/tests/et-legacy-comparison-contract.md](et-legacy-comparison-contract.md)
- [wiki/tests/expected-failure-and-error-keys.md](expected-failure-and-error-keys.md)
- [wiki/core/tests-and-fixtures.md](../core/tests-and-fixtures.md)
- [wiki/gems/flux-source/tests-and-fixtures.md](../gems/flux-source/tests-and-fixtures.md)
- [wiki/gems/con2prim/tests-and-fixtures.md](../gems/con2prim/tests-and-fixtures.md)
- [wiki/gems/eos/tests-and-fixtures.md](../gems/eos/tests-and-fixtures.md)
- [wiki/gems/induction/tests-and-fixtures.md](../gems/induction/tests-and-fixtures.md)
- [wiki/gems/neutrinos/tests-and-fixtures.md](../gems/neutrinos/tests-and-fixtures.md)
- [wiki/gems/reconstruction/tests-and-fixtures.md](../gems/reconstruction/tests-and-fixtures.md)
- [Unit_Tests/](../../Unit_Tests/)
- [Unit_Tests/data_gen/](../../Unit_Tests/data_gen/)
- [Unit_Tests/nrpyleakage_main.h](../../Unit_Tests/nrpyleakage_main.h)
- [Unit_Tests/sample_table/](../../Unit_Tests/sample_table/)
- [GRHayL/include/ghl_unit_tests.h](../../GRHayL/include/ghl_unit_tests.h)
- [GRHayL/include/make.code.defn](../../GRHayL/include/make.code.defn)
- [configure](../../configure)
- [.github/run_tests.sh](../../.github/run_tests.sh)
- [.github/workflows/](../../.github/workflows/)
