# Core Tests And Fixtures

This page routes Core-related tests, fixture generators, helper-only files, and
coverage gaps. It is a map for review impact; source, tests, and runner files
remain authoritative.

Use this page with [Test Map](../test-map.md),
[Con2Prim tests and fixtures](../gems/con2prim/tests-and-fixtures.md), and
[EOS tests and fixtures](../gems/eos/tests-and-fixtures.md). Do not store
binary metadata, hashes, `mtime`, fingerprints, or copied fixture contents here.

## Direct Core tests

### Core suite

- [Unit_Tests/unit_test_grhayl_core_test_suite.c](../../Unit_Tests/unit_test_grhayl_core_test_suite.c)
  is the direct Core suite. It checks simple EOS default handling, hybrid EOS
  setup, `ghl_set_prims_to_constant_atm` for simple and hybrid EOS, fixed-value
  4D vector-square and raise/lower helpers, metric determinant enforcement via
  `ghl_enforce_detgtij_and_initialize_ADM_metric`, and valid
  `ghl_get_con2prim_routine_name` keys.
- [Unit_Tests/unit_test_grhayl_core_test_suite.c](../../Unit_Tests/unit_test_grhayl_core_test_suite.c)
  opens `grhayl_core_test_suite_input.bin` and uses it as the metric fixture for
  determinant-enforced ADM metric initialization.
- [Unit_Tests/unit_test_grhayl_core_test_suite.c](../../Unit_Tests/unit_test_grhayl_core_test_suite.c)
  has a commented tabulated EOS atmosphere setup and a TODO before extending the
  atmosphere reset loop to tabulated EOS.
- Its simple/hybrid Atmosphere checks assert density, pressure, epsilon,
  entropy, and zero velocity. The `Y_e` and temperature assertions are guarded
  by `eos_type == 2`, which the loop never reaches, so this suite does not
  directly validate those two assignments for any EOS family. Its simple and
  hybrid EOS structs are not zero-initialized, and their initializers do not set
  those fields before the Atmosphere call reads them.

### Core fixture generator

- [Unit_Tests/data_gen/unit_test_data_grhayl_core_test_suite.c](../../Unit_Tests/data_gen/unit_test_data_grhayl_core_test_suite.c)
  generates random metric inputs through `ghl_randomize_metric` and writes the
  Core-suite metric fixture.
- Fixture naming drift is visible in repo-local files:
  [Unit_Tests/data_gen/unit_test_data_grhayl_core_test_suite.c](../../Unit_Tests/data_gen/unit_test_data_grhayl_core_test_suite.c)
  writes `grhayL_core_test_suite_input.bin`, while
  [Unit_Tests/unit_test_grhayl_core_test_suite.c](../../Unit_Tests/unit_test_grhayl_core_test_suite.c),
  [.github/run_tests.sh](../../.github/run_tests.sh), and
  [.github/workflows/](../../.github/workflows/) use/download
  `grhayl_core_test_suite_input.bin`.

### Direct helper checks visible in tests

- [Unit_Tests/unit_test_grhayl_core_test_suite.c](../../Unit_Tests/unit_test_grhayl_core_test_suite.c)
  directly checks `ghl_compute_vec2_from_vec4D` and
  `ghl_raise_lower_vector_4D` with fixed values.
- [Unit_Tests/unit_test_grhayl_core_test_suite.c](../../Unit_Tests/unit_test_grhayl_core_test_suite.c)
  exercises valid `ghl_get_con2prim_routine_name` enum keys; invalid-key
  behavior is covered in
  [Unit_Tests/unit_test_code_error.c](../../Unit_Tests/unit_test_code_error.c).
- [Unit_Tests/unit_test_code_error.c](../../Unit_Tests/unit_test_code_error.c)
  maps expected error codes for invalid EOS setup, invalid Fermi-Dirac key,
  singular `u0`, invalid Con2Prim keys, HDF5/table failures, invalid EOS/table
  state, and root-bracketing failures. Its key `33` checks invalid
  `ghl_get_con2prim_routine_name` behavior.
- [.github/run_tests.sh](../../.github/run_tests.sh) and
  [.github/workflows/](../../.github/workflows/) run
  `unit_test_code_error` as expected-failure keys, using those files only as
  fixture/test routing authority here.

## Indirect Core-dependent tests

### Stress-energy/conservatives

- [Unit_Tests/unit_test_compute_conservs_and_Tmunu.c](../../Unit_Tests/unit_test_compute_conservs_and_Tmunu.c)
  initializes Core metric, ADM auxiliaries, primitives, conservatives, and
  stress-energy structs; then validates `ghl_compute_conservs_and_Tmunu`,
  standalone `ghl_compute_conservs`, and standalone `ghl_compute_TDNmunu`
  against downloaded Con2Prim fixture data.
- [Unit_Tests/unit_test_compute_conservs_and_Tmunu.c](../../Unit_Tests/unit_test_compute_conservs_and_Tmunu.c)
  validates lower-index stress-energy fixture components through
  [Unit_Tests/pert_test_fail_stress_energy.c](../../Unit_Tests/pert_test_fail_stress_energy.c).
  It does not directly validate `ghl_compute_TUPmunu` against a fixture.
- [.github/run_tests.sh](../../.github/run_tests.sh) downloads
  `metric_Bfield_initial_data.bin` and
  `compute_conservs_and_Tmunu_{input,output,output_pert}.bin` from the
  Con2Prim fixture route before running this test.

### Velocity/u0

- [Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c](../../Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c)
  routes through Core metric initialization and ADM auxiliaries, then checks
  Con2Prim primitive limiting plus final `u0` values against fixtures.
- [Unit_Tests/unit_test_code_error.c](../../Unit_Tests/unit_test_code_error.c)
  directly exercises the Core `ghl_limit_v_and_compute_u0` singular-`u0` error
  path through key `4`.
- [.github/run_tests.sh](../../.github/run_tests.sh) downloads
  `metric_Bfield_initial_data.bin` and
  `enforce_primitive_limits_and_compute_u0_{input,output,output_pert}.bin` from
  the Con2Prim fixture route before running this test.

### ET legacy conservs

- [Unit_Tests/unit_test_ET_Legacy_conservs.c](../../Unit_Tests/unit_test_ET_Legacy_conservs.c)
  is an indirect Core regression route for metric initialization, ADM
  auxiliaries, primitive/conservative/stress-energy packing, primitive limiting,
  conservative/stress-energy computation, and `ghl_return_*` repacking in a
  legacy comparison path.
- [.github/run_tests.sh](../../.github/run_tests.sh) downloads
  `ET_Legacy_conservs_{input,output,output_pert}.bin` before running
  `unit_test_ET_Legacy_conservs`; workflows mirror that fixture route under
  [.github/workflows/](../../.github/workflows/).

### Broad Con2Prim/flux use of Core helpers

- [Unit_Tests/unit_test_con2prim_multi_method_hybrid.c](../../Unit_Tests/unit_test_con2prim_multi_method_hybrid.c),
  [Unit_Tests/unit_test_con2prim_tabulated.c](../../Unit_Tests/unit_test_con2prim_tabulated.c),
  [Unit_Tests/unit_test_hybrid_failure.c](../../Unit_Tests/unit_test_hybrid_failure.c),
  [Unit_Tests/unit_test_apply_conservative_limits.c](../../Unit_Tests/unit_test_apply_conservative_limits.c),
  and [Unit_Tests/unit_test_ET_Legacy_primitives.c](../../Unit_Tests/unit_test_ET_Legacy_primitives.c)
  repeatedly use Core initializers, metric setup, and ADM auxiliaries while
  testing Con2Prim behavior.
- [Unit_Tests/unit_test_hybrid_flux.c](../../Unit_Tests/unit_test_hybrid_flux.c)
  and [Unit_Tests/unit_test_tabulated_flux.c](../../Unit_Tests/unit_test_tabulated_flux.c)
  use Core metric, primitive initialization, and `ghl_limit_v_and_compute_u0`
  while testing flux behavior.
- [Unit_Tests/unit_test_ET_Legacy_flux_source.c](../../Unit_Tests/unit_test_ET_Legacy_flux_source.c)
  uses Core metric, extrinsic-curvature, and primitive initialization in a
  Flux_Source legacy comparison path.

## Helper-only files

### Installed test-support header

- [GRHayL/include/ghl_unit_tests.h](../../GRHayL/include/ghl_unit_tests.h) is
  named by [GRHayL/include/make.code.defn](../../GRHayL/include/make.code.defn),
  so it is in the configured install-header set.
- Some non-inline declarations are implemented by helper files under
  [Unit_Tests/](../../Unit_Tests/), including interpolation helpers,
  perturbation comparators, randomizers, and `get_table_quantity`. The declared
  binary read/write family and `ghl_initial_random_data` have no visible
  definitions. None of these test helpers appears in the Core/gem library
  source manifests, so install-header membership does not establish linkability
  from `libghl`; tests link only the helper sources selected by the test build.

### Perturbation helpers

- [Unit_Tests/pert_test_fail_conservatives.c](../../Unit_Tests/pert_test_fail_conservatives.c)
  compares conservative components against trusted and perturbed bars. It is a
  helper, not a standalone test.
- [Unit_Tests/pert_test_fail_primitives.c](../../Unit_Tests/pert_test_fail_primitives.c)
  compares primitive components, including entropy when enabled and tabulated
  `Y_e`/temperature when the EOS type is tabulated. It is a helper, not a
  standalone test.
- [Unit_Tests/pert_test_fail_stress_energy.c](../../Unit_Tests/pert_test_fail_stress_energy.c)
  compares named stress-energy components against trusted and perturbed bars.
  It is a helper, not a standalone test.

### Randomizers

- [Unit_Tests/randomize_metric.c](../../Unit_Tests/randomize_metric.c)
  builds random metric inputs with lapse fixed to `1.0` and shift fixed to
  zero. It feeds fixture generators, including the Core-suite generator.
- [Unit_Tests/randomize_primitives.c](../../Unit_Tests/randomize_primitives.c)
  builds random velocity and magnetic-field inputs around fixed Lorentz-factor
  and magnetic-pressure choices. It is helper code, not a standalone test.

## Missing/weak coverage

- Every pack/unpack round trip: weak direct coverage. The direct Core suite in
  [Unit_Tests/unit_test_grhayl_core_test_suite.c](../../Unit_Tests/unit_test_grhayl_core_test_suite.c)
  does not round-trip every `ghl_initialize_*`/`ghl_return_*` pair. The visible
  repack route is indirect in
  [Unit_Tests/unit_test_ET_Legacy_conservs.c](../../Unit_Tests/unit_test_ET_Legacy_conservs.c)
  for primitives, conservatives, and stress energy.
- Parameter initialization: callers use `ghl_initialize_params` throughout the
  suite, but no focused test in the listed Core ground truth asserts all copied
  fields, the `1e-8` default solver tolerance, the 30-iteration default, or all
  PPM defaults.
- EOS dispatch lifecycle: no focused test in the listed Core evidence calls a
  pointer before initialization, asserts all pointer targets after each EOS
  selection, or checks an invalid `ghl_eos_t` passed to the `void`
  `ghl_initialize_eos_functions` entry point.
- Simple EOS derived defaults: the Core suite asserts default `rho`/pressure
  floor and ceiling fields, but not the derived epsilon/entropy values produced
  after zero floors or other unchecked inputs such as `Gamma == 1`.
- Clamp helpers: no direct Core clamp test is visible in the listed ground
  truth. Keep `ghl_imin`, `ghl_imax`, `ghl_iclamp`, and `ghl_clamp` coverage as
  weak unless a direct test is added or found elsewhere in the repo.
- Standalone 3D metric helpers: direct fixed-value coverage in
  [Unit_Tests/unit_test_grhayl_core_test_suite.c](../../Unit_Tests/unit_test_grhayl_core_test_suite.c)
  is for the 4D helper paths. Standalone `ghl_raise_lower_vector_3D` and
  `ghl_compute_vec2_from_vec3D` fixture checks are not visible in the listed
  ground truth.
- Direct `ghl_compute_TUPmunu` fixture validation: not visible in
  [Unit_Tests/unit_test_compute_conservs_and_Tmunu.c](../../Unit_Tests/unit_test_compute_conservs_and_Tmunu.c),
  which validates combined conservative/stress-energy output and standalone
  `ghl_compute_TDNmunu` instead.
- Tabulated atmosphere reset branch: commented/TODO in
  [Unit_Tests/unit_test_grhayl_core_test_suite.c](../../Unit_Tests/unit_test_grhayl_core_test_suite.c)
  and routed as dependent EOS coverage in
  [EOS tests and fixtures](../gems/eos/tests-and-fixtures.md).
