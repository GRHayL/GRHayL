# Fixture Lifecycle And Harness Contract

This page routes Unit_Tests fixture lifecycle questions. Repo-local source files
remain ground truth for binary ordering, tolerance details, and executable
behavior.

Use this page for fixture family ownership, helper-only file classification,
generated-vs-downloaded data semantics, and drift routing. Do not use it as a
field-by-field binary specification.

## Source Authority

- Test inventory and owner routing: [test-map.md](../test-map.md).
- Drift routing: [contradictions.md](../contradictions.md).
- Runner download behavior: [run_tests.sh](../../.github/run_tests.sh) and
  workflow files under [.github/workflows](../../.github/workflows/).
- Harness tolerance helpers: [ghl_unit_tests.h](../../GRHayL/include/ghl_unit_tests.h),
  [pert_test_fail_conservatives.c](../../Unit_Tests/pert_test_fail_conservatives.c),
  [pert_test_fail_primitives.c](../../Unit_Tests/pert_test_fail_primitives.c),
  and [pert_test_fail_stress_energy.c](../../Unit_Tests/pert_test_fail_stress_energy.c).
- Fixture writers/readers: paired `Unit_Tests/unit_test_*.c` and
  `Unit_Tests/data_gen/unit_test_data_*.c` files.

## Harness Roles

- Trusted data: expected output read from downloaded TestData fixtures or from
  generated local output files. Tests compare computed values against this data.
- Computed data: values produced by the test run from GRHayL routines and test
  helper calls.
- Perturbed data: comparison bar data, usually `*_output_pert.bin` or
  `*_perturbed.bin`; exact perturbation magnitude and exact tolerance behavior
  belong in source.

Most replay tests use `ghl_pert_test_fail` or wrapper helpers around it. The
exact relative-error cutoff, absolute cutoff, custom cutoffs, and per-quantity
exceptions route to `GRHayL/include/ghl_unit_tests.h` plus the Unit_Tests helper
source above. Wiki pages should describe roles, not duplicate tolerance logic.

## Helper-Only Files

These files are helper-only evidence. They are compiled into tests or data
generators but are not standalone unit tests:

| Helper-only file | Role | Route |
| --- | --- | --- |
| `pert_test_fail_conservatives.c` | Conservative quantity comparison wrapper around `ghl_pert_test_fail`. | [Core tests and fixtures](../core/tests-and-fixtures.md), [Con2Prim tests and fixtures](../gems/con2prim/tests-and-fixtures.md) |
| `pert_test_fail_primitives.c` | Primitive quantity comparison wrapper, including EOS-dependent fields and cutoffs. | [Con2Prim tests and fixtures](../gems/con2prim/tests-and-fixtures.md), [EOS tests and fixtures](../gems/eos/tests-and-fixtures.md) |
| `pert_test_fail_stress_energy.c` | Stress-energy tensor comparison wrapper. | [Core tests and fixtures](../core/tests-and-fixtures.md) |
| `randomize_metric.c` | Random metric setup used by generators/tests. | [Core tests and fixtures](../core/tests-and-fixtures.md) |
| `randomize_primitives.c` | Random primitive setup used by generators/tests. | [Con2Prim tests and fixtures](../gems/con2prim/tests-and-fixtures.md) |
| `test_compute_h_and_cs2.c` | EOS enthalpy/sound-speed check helper. | [EOS tests and fixtures](../gems/eos/tests-and-fixtures.md) |
| `tabulated_eos_unit_test_helpers.c` | Analytic table quantity helper for tabulated EOS checks. | [EOS tests and fixtures](../gems/eos/tests-and-fixtures.md) |
| `compute_A_flux_with_B.c` | Induction vector-potential HLL helper using `B`. | [Induction tests and fixtures](../gems/induction/tests-and-fixtures.md) |
| `compute_A_flux_with_Btilde.c` | Induction vector-potential HLL helper using `Btilde`. | [Induction tests and fixtures](../gems/induction/tests-and-fixtures.md) |
| `compute_ccc_ADM.c` | Cell-centered ADM interpolation helper. | [Induction tests and fixtures](../gems/induction/tests-and-fixtures.md) |
| `compute_ccc_BSSN.c` | Cell-centered BSSN interpolation helper. | [Induction tests and fixtures](../gems/induction/tests-and-fixtures.md) |
| `compute_vvv_ADM.c` | Vertex-centered ADM interpolation helper. | [Induction tests and fixtures](../gems/induction/tests-and-fixtures.md) |

`nrpyleakage_main.h` is also helper-only: it owns shared NRPyLeakage main/key
handling and routes through [Neutrinos tests and fixtures](../gems/neutrinos/tests-and-fixtures.md).

## Fixture Source Classes

| Source class | Meaning | Evidence |
| --- | --- | --- |
| Local data generators | `Unit_Tests/data_gen/unit_test_data_*.c` programs write local fixture files when built and run through data-generation paths. | [Unit_Tests/data_gen](../../Unit_Tests/data_gen/) |
| Generation modes in tests | Some test binaries generate their own `*_unperturbed.bin` and `*_perturbed.bin` families when invoked with a generation key, then replay when invoked with test key. | NRPyLeakage tests and tabulated Con2Prim tests |
| Downloaded TestData fixtures | Normal CI/script runs download binary fixtures from `GRHayL/TestData` into the repo root before running tests. | [run_tests.sh](../../.github/run_tests.sh), workflows |
| Downloaded EOS tables | Scripted runs download large HDF5 tables from `stellarcollapse.org/EOS`; these are test inputs, not generated binary bars. | [run_tests.sh](../../.github/run_tests.sh) |
| Checked-in sample tables | Reduced/analytic sample table assets live under `Unit_Tests/sample_table/`; they are input assets, not replay output bars. | [EOS tests and fixtures](../gems/eos/tests-and-fixtures.md) |

Downloaded `TestData` fixtures and locally generated fixtures may share names.
Normal runner behavior uses downloads; generator presence does not mean CI
regenerates trusted data.

## Fixture Families By Owner

| Owner | Fixture families | Notes |
| --- | --- | --- |
| Core/chalice | `grhayl_core_test_suite_input.bin` | Generator spelling drift with `grhayL_core_test_suite_input.bin` routes to [contradictions.md](../contradictions.md). |
| Con2Prim/Core shared | `metric_Bfield_initial_data.bin`, `apply_conservative_limits_*`, `con2prim_multi_method_hybrid_*`, `enforce_primitive_limits_and_compute_u0_*`, `compute_conservs_and_Tmunu_*` | `metric_Bfield_initial_data.bin` is shared by multiple tests; changing it has broad blast radius across limit, recovery, `u0`, conservative, and stress-energy checks. |
| Con2Prim tabulated | `con2prim_tabulated_*_unperturbed.bin`, `con2prim_tabulated_*_perturbed.bin` | HDF5-table-backed replay/generation modes route through Con2Prim/EOS pages. |
| EOS | `simple_table.h5`, checked-in sample-table assets | Table lifecycle routes through EOS test docs. |
| Flux_Source | `hybrid_flux_*`, `tabulated_flux_*` | Tabulated flux also depends on downloaded LS220 EOS table. |
| Induction | `induction_interpolation_*`, `HLL_flux_*` | HLL here is vector-potential induction evidence, not Flux_Source HLLE evidence. |
| Neutrinos | `nrpyleakage_optically_thin_gas_*`, `nrpyleakage_constant_density_sphere_*`, `nrpyleakage_luminosities_*` | Families use `*_unperturbed.bin`/`*_perturbed.bin` and SLy4 table input. |
| Reconstruction | `PLM_reconstruction_*`, `WENOZ_reconstruction_*`, `ET_Legacy_reconstruction_*` | ET Legacy reconstruction exercises PPM legacy comparison paths. |
| ET Legacy | `ET_Legacy_conservs_*`, `ET_Legacy_primitives_*`, `ET_Legacy_flux_source_*`, `ET_Legacy_HLL_flux_*`, `ET_Legacy_induction_gauge_rhs_*`, `ET_Legacy_reconstruction_*` | See [ET legacy comparison contract](et-legacy-comparison-contract.md). |

Binary byte layout belongs in paired test/generator source. Wiki may mention
`fread`/`fwrite` only to route readers back to source; it must not publish
field-by-field layout tables.

## Drift Routing

- `Test_Data` vs `TestData` fixture repository spelling routes to
  [contradictions.md](../contradictions.md).
- Fixed perturbation wording vs test-specific perturbation behavior routes to
  [contradictions.md](../contradictions.md).
- `grhayL_core_test_suite_input.bin` vs `grhayl_core_test_suite_input.bin`
  fixture spelling routes to [contradictions.md](../contradictions.md).
