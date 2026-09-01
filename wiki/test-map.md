# Test Map

Repo ground truth: `Unit_Tests/`, `.github/run_tests.sh`, `.github/workflows/`,
`configure`, `README.md`, and `docs/raw/mainpage.md`.

Evidence labels are strict: `configure` selects targets; `make tests` and
`make datagen` compile/link; only an exact invocation establishes execution.
Workflow commands are workflow-only evidence, not historical pass results.
Default configuration selects 30 unit-test binaries; `.github/run_tests.sh`
directly invokes 27, omitting WENOZ reconstruction, Con2Prim debug, and the
CompOSE integration test. WENOZ and CompOSE are workflow-selected; no normal
invocation for the debug binary is visible.

Core/chalice test selection, fixture naming, helper-only files, and weak
coverage notes route through [Core tests and fixtures](core/tests-and-fixtures.md).
Con2Prim-specific test selection, fixture names, HDF5 gates, and debug-test
caveats route through [Con2Prim tests and fixtures](gems/con2prim/tests-and-fixtures.md).
EOS-specific direct tests, helper-only files, sample table roles, generated
`simple_table.h5`, and no-HDF5 error coverage route through
[EOS tests and fixtures](gems/eos/tests-and-fixtures.md).
Neutrinos-specific CLI keys, fixture pairs, SLy4 table setup, HDF5 gates, and
CI downloads route through [Neutrinos tests and fixtures](gems/neutrinos/tests-and-fixtures.md).
Reconstruction-specific PLM, WENOZ, ET Legacy fixtures, generators, CI nuance,
and PPM coverage caveats route through [Reconstruction tests and fixtures](gems/reconstruction/tests-and-fixtures.md).
Induction-specific fixture/test details route through [Induction tests and fixtures](gems/induction/tests-and-fixtures.md);
targeted build, run, and CI job guidance routes through [Induction verification workflows](gems/induction/verification-workflows.md).
Flux_Source-specific hybrid/tabulated HLLE fixtures, ET Legacy flux/source
replay, characteristic-speed fixture evidence, and the Induction HLL exclusion
route through [Flux_Source tests and fixtures](gems/flux-source/tests-and-fixtures.md).

Cross-cutting Unit_Tests routes live under [Unit_Tests Hub](tests/index.md):
[runner and generated artifacts](tests/runner-and-generated-artifacts.md),
[fixture lifecycle and harness contract](tests/fixture-lifecycle-and-harness-contract.md),
[HDF5 sample tables](tests/hdf5-sample-tables.md),
[ET Legacy comparison contract](tests/et-legacy-comparison-contract.md),
[expected-failure and error keys](tests/expected-failure-and-error-keys.md),
and [unit-test coverage and gap matrix](tests/unit-test-coverage-and-gap-matrix.md).
This page remains the top-level inventory; use the detailed pages for fixture
lifecycle, run modes, tables, and coverage caveats.

GRHayLib/Cactus thorn coverage caveat: ET_Legacy tests are related legacy
comparisons for upstream GRHayL behavior, not direct GRHayLib thorn coverage.
Implementation-only changes under `implementations/GRHayLib/` need manual
downstream verification, including Cactus or Einstein Toolkit checks when
available. Route that checklist through
[GRHayLib verification and drift](implementations/grhaylib/verification-and-drift.md).

## Unit Tests

| Test | Likely source area | Behavior covered | Fixture dependencies | Notes/gaps |
| --- | --- | --- | --- | --- |
| `Unit_Tests/unit_test_ET_Legacy_HLL_flux.c` | `GRHayL/Induction/`, legacy comparison | Vector-potential HLL flux against ET legacy outputs; see [Induction HLL flux contract](gems/induction/hll-flux-contract.md). | `ET_Legacy_HLL_flux_{input,output,output_pert}.bin` from repo-visible `GRHayL/TestData` path; route fixture details through [Induction tests and fixtures](gems/induction/tests-and-fixtures.md). | Data generator writes input and perturbed input only; trusted output comes from external legacy path per docs. |
| `Unit_Tests/unit_test_ET_Legacy_conservs.c` | `GRHayL/Con2Prim/`, `GRHayL/GRHayL_Core/` | Primitive-to-conservative legacy comparison; see [limits and conversions](gems/con2prim/limits-and-conversions.md) and [Core tests and fixtures](core/tests-and-fixtures.md). | `ET_Legacy_conservs_{input,output,output_pert}.bin`. | Perturbation magnitude differs by generator; see `contradictions.md`. |
| `Unit_Tests/unit_test_ET_Legacy_flux_source.c` | `GRHayL/Flux_Source/` | Flux/source RHS legacy comparison; see [Flux_Source tests and fixtures](gems/flux-source/tests-and-fixtures.md). | `ET_Legacy_flux_source_{input,output,output_pert}.bin`. | Legacy output source is external by docs. |
| `Unit_Tests/unit_test_ET_Legacy_induction_gauge_rhs.c` | `GRHayL/Induction/` | BSSN interpolation plus vector-potential and gauge RHS comparison; see [Induction gauge RHS contract](gems/induction/gauge-rhs-contract.md). | `ET_Legacy_induction_gauge_rhs_{input,output,output_pert}.bin`; route fixture details through [Induction tests and fixtures](gems/induction/tests-and-fixtures.md). | Uses helper interpolation behavior internally. |
| `Unit_Tests/unit_test_ET_Legacy_primitives.c` | `GRHayL/Con2Prim/` | Conservative-to-primitive legacy comparison; see [recovery flow](gems/con2prim/recovery-flow.md). | `ET_Legacy_primitives_{input,output,output_pert}.bin`. | Generator uses a larger perturbation than some other ET legacy generators. |
| `Unit_Tests/unit_test_ET_Legacy_reconstruction.c` | `GRHayL/Reconstruction/` | Legacy reconstruction comparison; see [Reconstruction tests and fixtures](gems/reconstruction/tests-and-fixtures.md). | `ET_Legacy_reconstruction_{input,output,output_pert}.bin`. | Generator writes input and perturbed input only; replay does not read perturbed input. ET Legacy exercises PPM paths. |
| `Unit_Tests/unit_test_HLL_flux.c` | `GRHayL/Induction/` | HLL vector-potential flux with `B` and `Btilde`; see [Induction HLL flux contract](gems/induction/hll-flux-contract.md). | `HLL_flux_input.bin`, `HLL_flux_with_B_output*.bin`, `HLL_flux_with_Btilde_output*.bin`; route fixture details through [Induction tests and fixtures](gems/induction/tests-and-fixtures.md). | Shared helper functions compute expected variants during generation. |
| `Unit_Tests/unit_test_PLM_reconstruction.c` | `GRHayL/Reconstruction/PLM/` | Minmod, MC, and superbee PLM right/left face reconstruction; see [PLM limiters](gems/reconstruction/plm-limiters.md) and [Reconstruction tests and fixtures](gems/reconstruction/tests-and-fixtures.md). | `PLM_reconstruction_{input,output,output_pert}.bin`. | Does not cover PPM. |
| `Unit_Tests/unit_test_WENOZ_reconstruction.c` | `GRHayL/Reconstruction/WENOZ/` | WENOZ reconstruction right/left face values; see [WENOZ contract](gems/reconstruction/wenoz-contract.md) and [Reconstruction tests and fixtures](gems/reconstruction/tests-and-fixtures.md). | `WENOZ_reconstruction_{input,output,output_pert}.bin`. | No dedicated `unit_test_PPM_reconstruction.c` visible in repo file list. |
| `Unit_Tests/unit_test_apply_conservative_limits.c` | `GRHayL/Con2Prim/` | Conservative floor/limit application; see [limits and conversions](gems/con2prim/limits-and-conversions.md). | `metric_Bfield_initial_data.bin`, `apply_conservative_limits_{input,output,output_pert}.bin`. | Uses shared metric/B-field fixture. |
| `Unit_Tests/unit_test_code_error.c` | error handling across modules | Expected error-code paths for EOS, C2P, velocity limiting, tabulated EOS, NN loader/init, and NRPyLeakage; see [expected-failure and error keys](tests/expected-failure-and-error-keys.md), [Core tests and fixtures](core/tests-and-fixtures.md), and [EOS tests and fixtures](gems/eos/tests-and-fixtures.md). | HDF5-enabled runs need `SLy4_3335_rho391_temp163_ye66.h5` for many tabulated setup/interpolation/root keys; table-read failure keys use temporary or missing `test.h5`; NN key `85` uses the same SLy4 table without an embedded `grhayl_nn_c2p` group. | `run_tests.sh` expects process failure for keys `0` through `85`; HDF5-only keys are skipped in no-HDF5 builds. |
| `Unit_Tests/unit_test_compute_conservs_and_Tmunu.c` | `GRHayL/Con2Prim/`, `GRHayL/GRHayL_Core/` | Conservative variables and stress-energy tensor, combined and standalone calls; see [Core tests and fixtures](core/tests-and-fixtures.md) and [limits and conversions](gems/con2prim/limits-and-conversions.md). | `metric_Bfield_initial_data.bin`, `compute_conservs_and_Tmunu_{input,output,output_pert}.bin`. | Compares against perturbation bars. |
| `Unit_Tests/unit_test_con2prim_debug.c` | `GRHayL/Con2Prim/Tabulated/` | Single hard-coded tabulated C2P debug case. | CLI EOS table path. | Not run by `run_tests.sh`; excluded from no-HDF5 configured tests. See [tests and fixtures](gems/con2prim/tests-and-fixtures.md). |
| `Unit_Tests/unit_test_con2prim_multi_method_hybrid.c` | `GRHayL/Con2Prim/Hybrid/` | Hybrid selected-method C2P recovery for supported hybrid solver keys; see [solver matrix](gems/con2prim/solver-matrix.md) and [recovery flow](gems/con2prim/recovery-flow.md). | `metric_Bfield_initial_data.bin`, `con2prim_multi_method_hybrid_{input,output,output_pert}.bin`. | Selector gets undensitized `cons_undens`, but separate `ghl_guess_primitives` check receives densitized `cons`; that call does not validate helper's undensitized contract. |
| `Unit_Tests/unit_test_con2prim_tabulated.c` | `GRHayL/Con2Prim/Tabulated/`, `GRHayL/EOS/Tabulated/` | Tabulated multi-method C2P for Palenzuela1D, Newman1D, and Noble2D style routines; with nonzero key, replays with NN primitive guesses disabled and enabled; see [solver matrix](gems/con2prim/solver-matrix.md). | CLI EOS table path plus generated/downloaded `con2prim_tabulated_*_{unperturbed,perturbed}.bin`; runner/workflows run `pyghl append SLy4_3335_rho391_temp163_ye66.h5` before key `1`. | Generates data when invoked with generation mode in source; `run_tests.sh` downloads fixtures, records only the visible `pyghl append` setup command, and runs key `1`. See [tests and fixtures](gems/con2prim/tests-and-fixtures.md). |
| `Unit_Tests/unit_test_c2p_nn_guess.c` | `GRHayL/Con2Prim/Tabulated/neural_network_guess/`, `GRHayL/EOS/Tabulated/` | Direct NN primitive-guess model validation, inference fallback, root/embedded/legacy HDF5 loads, malformed HDF5 paths, and failed-load preservation. | No external binary fixture. HDF5-enabled runs create temporary `/tmp/unit_test_c2p_nn_*.h5` files. | Run by `run_tests.sh` and workflow `c2p-failure` matrix entry `c2p_nn_guess`; HDF5 loader checks are compiled out under `GHL_DISABLE_HDF5`. |
| `Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c` | `GRHayL/Con2Prim/`, `GRHayL/GRHayL_Core/` | Primitive floors/ceilings and `u0` computation; see [Core tests and fixtures](core/tests-and-fixtures.md) and [limits and conversions](gems/con2prim/limits-and-conversions.md). | `metric_Bfield_initial_data.bin`, `enforce_primitive_limits_and_compute_u0_{input,output,output_pert}.bin`. | Shared metric/B-field fixture. |
| `Unit_Tests/unit_test_grhayl_core_test_suite.c` | `GRHayL/GRHayL_Core/` | Simple EOS default bounds, hybrid EOS setup, constant Atmosphere reset for simple/hybrid EOS, fixed 4D metric helpers, determinant-enforced ADM metric initialization, and valid Con2Prim routine-name keys; see [Core tests and fixtures](core/tests-and-fixtures.md). | `grhayl_core_test_suite_input.bin`. | Data generator string uses `grhayL_core_test_suite_input.bin`; script downloads lowercase `grhayl_core_test_suite_input.bin`. Tabulated Atmosphere reset branch is commented/TODO; see [EOS tests and fixtures](gems/eos/tests-and-fixtures.md). |
| `Unit_Tests/unit_test_hybrid_failure.c` | `GRHayL/Con2Prim/Hybrid/` | Expected hybrid C2P failure codes and backup routine behavior; see [recovery flow](gems/con2prim/recovery-flow.md). | None visible; hard-coded cases. | No external fixture. |
| `Unit_Tests/unit_test_hybrid_flux.c` | `GRHayL/Flux_Source/hybrid*` | HLLE fluxes for hybrid EOS, with entropy and three directions; see [Flux_Source tests and fixtures](gems/flux-source/tests-and-fixtures.md). | `hybrid_flux_{input,output,output_pert}.bin`. | Covers generated/derived flux functions through public entry points. |
| `Unit_Tests/unit_test_induction_ccc_ADM.c` | `GRHayL/Induction/Interpolators/` | Cell-centered ADM interpolation; see [Induction interpolation and staggering contract](gems/induction/interpolation-and-staggering-contract.md). | Fixture family `induction_interpolation_*`; details in [Induction tests and fixtures](gems/induction/tests-and-fixtures.md). | Uses helper implementation in `Unit_Tests/compute_ccc_ADM.c`. |
| `Unit_Tests/unit_test_induction_ccc_BSSN.c` | `GRHayL/Induction/Interpolators/` | Cell-centered BSSN interpolation; see [Induction interpolation and staggering contract](gems/induction/interpolation-and-staggering-contract.md). | `induction_interpolation_input.bin`, `induction_interpolation_BSSN_input.bin`, `induction_interpolation_ccc_BSSN_output*.bin`; route fixture details through [Induction tests and fixtures](gems/induction/tests-and-fixtures.md). | Uses helper implementation in `Unit_Tests/compute_ccc_BSSN.c`. |
| `Unit_Tests/unit_test_induction_vvv_ADM.c` | `GRHayL/Induction/Interpolators/` | Vertex-centered ADM interpolation; see [Induction interpolation and staggering contract](gems/induction/interpolation-and-staggering-contract.md). | `induction_interpolation_input.bin`, `induction_interpolation_ADM_input.bin`, `induction_interpolation_vvv_ADM_output*.bin`; route fixture details through [Induction tests and fixtures](gems/induction/tests-and-fixtures.md). | Uses helper implementation in `Unit_Tests/compute_vvv_ADM.c`. |
| `Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c` | `GRHayL/Neutrinos/NRPyLeakage/` | Constant-density sphere opacities and optical-depth iteration; see [tests and fixtures](gems/neutrinos/tests-and-fixtures.md). | CLI EOS table path plus `nrpyleakage_constant_density_sphere_{unperturbed,perturbed}.bin`. | Comparison result is discarded; trusted/perturbed and neighbor-direction arguments are reversed. Execution is not numerical-pass evidence. |
| `Unit_Tests/unit_test_nrpyleakage_luminosities.c` | `GRHayL/Neutrinos/NRPyLeakage/` | Fermi-Dirac branch checks and neutrino luminosity replay; see [tests and fixtures](gems/neutrinos/tests-and-fixtures.md). | CLI EOS table path plus `nrpyleakage_luminosities_{unperturbed,perturbed}.bin`. | Comparison result is discarded, so numerical mismatch cannot fail executable. |
| `Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c` | `GRHayL/Neutrinos/NRPyLeakage/` | Optically thin gas leakage source evolution; see [tests and fixtures](gems/neutrinos/tests-and-fixtures.md). | CLI EOS table path plus `nrpyleakage_optically_thin_gas_{unperturbed,perturbed}.bin`. | Comparison result is discarded, so numerical mismatch cannot fail executable. |
| `Unit_Tests/unit_test_piecewise_polytrope.c` | `GRHayL/EOS/Hybrid/` | Piecewise-polytrope `K_ppoly` and `eps_integ_const` setup; see [EOS tests and fixtures](gems/eos/tests-and-fixtures.md). | None visible. | No external fixture. |
| `Unit_Tests/unit_test_tabulated_eos.c` | `GRHayL/EOS/Tabulated/` | HDF5 table read, analytic table quantity checks, tabulated interpolation routines, bounds, `ghl_compute_h_and_cs2`, and beta-equilibrium helpers; see [EOS tests and fixtures](gems/eos/tests-and-fixtures.md). | CLI table path, normally `simple_table.h5` from CI or local sample generator. | HDF5-only. |
| `Unit_Tests/unit_test_tabulated_eos_compose.c` | `tools/compose/`, unchanged tabulated EOS, Con2Prim, Flux_Source, and NRPyLeakage runtimes | Reads every node/all 19 fields independently; checks loader units/order, enthalpy, relativistic sound speed, midpoint interpolation, `eps/P/S/h` inverses from distinct valid initial guesses, and six-output order/ranges; requires no-fallback Palenzuela recovery; compares characteristic speed, HLLE/entropy fluxes, and source terms with analytic goldens; compares all eight combined-leakage outputs with fixed regression goldens for the two qualified table dimensions; and checks cleanup. Python failure and 100% branch coverage lives under `Unit_Tests/compose/`. | CLI regularized StellarCollapse table path; focused Ubuntu GCC CI builds an asymmetric synthetic table, while full table 141 is external manual qualification data. | HDF5-only; converter output must pass with runtime sound-speed cleaning disabled. |
| `Unit_Tests/unit_test_tabulated_flux.c` | `GRHayL/Flux_Source/tabulated*`, `GRHayL/EOS/Tabulated/` | HLLE fluxes for tabulated EOS, with entropy and three directions; see [Flux_Source tests and fixtures](gems/flux-source/tests-and-fixtures.md). | `LS220_234r_136t_50y_analmu_20091212_SVNr26.h5`, `tabulated_flux_{input,output,output_pert}.bin`. | HDF5-only; EOS table downloaded by `run_tests.sh`. |

## Data Generators

| Data generator | Likely source area | Behavior covered | Fixture dependencies | Notes/gaps |
| --- | --- | --- | --- | --- |
| `Unit_Tests/data_gen/unit_test_data_ET_Legacy_HLL_flux.c` | `GRHayL/Induction/` | Generates ET legacy HLL flux input and perturbed input; see [Induction tests and fixtures](gems/induction/tests-and-fixtures.md). | None visible beyond generated files. | Trusted output is external legacy data. |
| `Unit_Tests/data_gen/unit_test_data_ET_Legacy_conservs.c` | `GRHayL/Con2Prim/`, core | Generates primitive/conservative ET legacy input and perturbation; see [Core tests and fixtures](core/tests-and-fixtures.md). | None visible beyond generated files. | Uses random relative perturbations. |
| `Unit_Tests/data_gen/unit_test_data_ET_Legacy_flux_source.c` | `GRHayL/Flux_Source/` | Generates ET legacy flux/source input and perturbation; see [Flux_Source tests and fixtures](gems/flux-source/tests-and-fixtures.md). | None visible beyond generated files. | Trusted output is external legacy data. |
| `Unit_Tests/data_gen/unit_test_data_ET_Legacy_induction_gauge_rhs.c` | `GRHayL/Induction/` | Generates ET legacy induction gauge RHS input and perturbation; see [Induction tests and fixtures](gems/induction/tests-and-fixtures.md). | None visible beyond generated files. | Trusted output is external legacy data. |
| `Unit_Tests/data_gen/unit_test_data_ET_Legacy_primitives.c` | `GRHayL/Con2Prim/` | Generates ET legacy primitive recovery input and perturbation. | None visible beyond generated files. | Perturbation magnitude differs from README fixed wording. |
| `Unit_Tests/data_gen/unit_test_data_ET_Legacy_reconstruction.c` | `GRHayL/Reconstruction/` | Generates ET legacy reconstruction input and perturbation; see [Reconstruction tests and fixtures](gems/reconstruction/tests-and-fixtures.md). | None visible beyond generated files. | Trusted output is external legacy data. |
| `Unit_Tests/data_gen/unit_test_data_HLL_flux.c` | `GRHayL/Induction/` | Generates HLL flux inputs and output bars for `B` and `Btilde` variants; see [Induction HLL flux contract](gems/induction/hll-flux-contract.md) and [tests and fixtures](gems/induction/tests-and-fixtures.md). | None visible beyond generated files. | Local generator covers expected output. |
| `Unit_Tests/data_gen/unit_test_data_PLM_reconstruction.c` | `GRHayL/Reconstruction/PLM/` | Intended PLM input/output bars for minmod, MC, and superbee; see [Reconstruction tests and fixtures](gems/reconstruction/tests-and-fixtures.md). | None visible beyond generated files. | Boundary indexing can exceed arrays; do not call generator safe until fixed. No PPM generator visible. |
| `Unit_Tests/data_gen/unit_test_data_WENOZ_reconstruction.c` | `GRHayL/Reconstruction/WENOZ/` | Intended WENOZ input/output bars; see [Reconstruction tests and fixtures](gems/reconstruction/tests-and-fixtures.md). | None visible beyond generated files. | Boundary indexing can exceed arrays; do not call generator safe until fixed. |
| `Unit_Tests/data_gen/unit_test_data_con2prim_multi_method_hybrid.c` | `GRHayL/Con2Prim/Hybrid/` | Generates hybrid C2P plus shared metric/B-field and limit fixtures; see [tests and fixtures](gems/con2prim/tests-and-fixtures.md). | None visible beyond generated files. | Also writes inputs used by conservative limit and primitive limit tests. |
| `Unit_Tests/data_gen/unit_test_data_grhayl_core_test_suite.c` | `GRHayL/GRHayL_Core/` | Generates core suite input fixture; see [Core tests and fixtures](core/tests-and-fixtures.md). | None visible beyond generated files. | Output spelling in source appears as `grhayL_core_test_suite_input.bin`; downloaded fixture path uses `grhayl_core_test_suite_input.bin`. |
| `Unit_Tests/data_gen/unit_test_data_hybrid_flux.c` | `GRHayL/Flux_Source/hybrid*` | Generates hybrid HLLE flux inputs and output bars; see [Flux_Source tests and fixtures](gems/flux-source/tests-and-fixtures.md). | None visible beyond generated files. | Covers entropy and directional combinations. |
| `Unit_Tests/data_gen/unit_test_data_induction_interpolation.c` | `GRHayL/Induction/Interpolators/` | Generates interpolation inputs and output bars for ADM/BSSN cell/vertex cases; see [Induction interpolation and staggering contract](gems/induction/interpolation-and-staggering-contract.md) and [tests and fixtures](gems/induction/tests-and-fixtures.md). | None visible beyond generated files. | Feeds three induction interpolation tests. |
| `Unit_Tests/data_gen/unit_test_data_tabulated_flux.c` | `GRHayL/Flux_Source/tabulated*` | Generates tabulated HLLE flux inputs and output bars; see [Flux_Source tests and fixtures](gems/flux-source/tests-and-fixtures.md). | Requires `LS220_234r_136t_50y_analmu_20091212_SVNr26.h5`. | HDF5-only. |

## Sample EOS Table Files

| File | Role | Notes/gaps |
| --- | --- | --- |
| `Unit_Tests/sample_table/Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817_simple.h5` | Checked-in reduced HDF5 EOS table sample; see [EOS tests and fixtures](gems/eos/tests-and-fixtures.md). | Source provenance is described in `table_info.txt`; avoid relying on local absolute paths recorded there. |
| `Unit_Tests/sample_table/generate_simple_table.py` | Generates `simple_table.h5` with analytic EOS quantities for tabulated EOS tests; see [EOS tests and fixtures](gems/eos/tests-and-fixtures.md). | Uses `numpy` and `h5py`; output name differs from checked-in Hempel sample name. |
| `Unit_Tests/sample_table/table_info.txt` | Human-readable dimensions/ranges for input and reduced output table; see [EOS tests and fixtures](gems/eos/tests-and-fixtures.md). | Contains local absolute paths from the original table generation environment; do not copy those paths into docs. |

## Unit-Test Helpers

| Helper | Likely source area | Behavior covered | Fixture dependencies | Notes/gaps |
| --- | --- | --- | --- | --- |
| `Unit_Tests/compute_A_flux_with_B.c` | `GRHayL/Induction/` | Test helper for vector-potential HLL flux using `B`; see [Induction HLL flux contract](gems/induction/hll-flux-contract.md). | Used by `unit_test_HLL_flux.c` and generator. | Helper, not standalone test; route through [Induction tests and fixtures](gems/induction/tests-and-fixtures.md). |
| `Unit_Tests/compute_A_flux_with_Btilde.c` | `GRHayL/Induction/` | Test helper for vector-potential HLL flux using `Btilde`; see [Induction HLL flux contract](gems/induction/hll-flux-contract.md). | Used by `unit_test_HLL_flux.c` and generator. | Helper, not standalone test; route through [Induction tests and fixtures](gems/induction/tests-and-fixtures.md). |
| `Unit_Tests/compute_ccc_ADM.c` | `GRHayL/Induction/Interpolators/` | Cell-centered ADM interpolation helper; see [Induction interpolation and staggering contract](gems/induction/interpolation-and-staggering-contract.md). | Used by induction interpolation tests. | Helper, not standalone test; route through [Induction tests and fixtures](gems/induction/tests-and-fixtures.md). |
| `Unit_Tests/compute_ccc_BSSN.c` | `GRHayL/Induction/Interpolators/` | Cell-centered BSSN interpolation helper; see [Induction interpolation and staggering contract](gems/induction/interpolation-and-staggering-contract.md). | Used by induction interpolation tests. | Helper, not standalone test; route through [Induction tests and fixtures](gems/induction/tests-and-fixtures.md). |
| `Unit_Tests/compute_vvv_ADM.c` | `GRHayL/Induction/Interpolators/` | Vertex-centered ADM interpolation helper; see [Induction interpolation and staggering contract](gems/induction/interpolation-and-staggering-contract.md). | Used by induction interpolation tests. | Helper, not standalone test; route through [Induction tests and fixtures](gems/induction/tests-and-fixtures.md). |
| `Unit_Tests/nrpyleakage_main.h` | `GRHayL/Neutrinos/NRPyLeakage/` | Shared main/argument handling for NRPyLeakage tests; see [Neutrinos tests and fixtures](gems/neutrinos/tests-and-fixtures.md). | CLI EOS table path. | Header helper owns key `0` generation mode and key `1` unit-test mode. |
| `Unit_Tests/pert_test_fail_conservatives.c` | test harness | Conservative quantity perturbation-bar comparison. | Called by multiple tests. | Helper, not standalone test. |
| `Unit_Tests/pert_test_fail_primitives.c` | test harness | Primitive quantity perturbation-bar comparison, including EOS-dependent fields. | Called by multiple tests. | Helper, not standalone test. |
| `Unit_Tests/pert_test_fail_stress_energy.c` | test harness | Stress-energy perturbation-bar comparison. | Called by stress-energy tests. | Helper, not standalone test. |
| `Unit_Tests/randomize_metric.c` | test harness | Random metric generation. | Used by data generators/tests. | Helper, not standalone test. |
| `Unit_Tests/randomize_primitives.c` | test harness | Random primitive generation. | Used by data generators/tests. | Helper, not standalone test. |
| `Unit_Tests/tabulated_eos_unit_test_helpers.c` | `GRHayL/EOS/Tabulated/` | Analytic helper functions for tabulated EOS validation. | Used by tabulated EOS tests. | HDF5-adjacent helper. |
| `Unit_Tests/test_compute_h_and_cs2.c` | `GRHayL/EOS/` | Helper for enthalpy and sound-speed checks. | Used by EOS/flux tests. | Helper, not standalone test. |

## Repo-Visible External Fixture References

`README.md` names `GRHayL/Test_Data`, while `docs/raw/mainpage.md`,
`.github/run_tests.sh`, and workflows use `GRHayL/TestData`. The scripts also
reference EOS tables under `stellarcollapse.org/EOS`. This page records only
references visible in repo files; no external lookup was used here.

Installed `GRHayL/include/ghl_unit_tests.h` remains test-only/unresolved
surface: inline helpers compile in callers, but non-inline declarations are not
generally defined by `libghl`. See [Public API Map](public-api-map.md).
