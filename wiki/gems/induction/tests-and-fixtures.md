# Induction Tests And Fixtures

## Purpose

This page routes Induction test binaries, helper-only files, data generators,
fixture ownership, and known coverage gaps. Repo-local tests, generators,
[Test Map](../../test-map.md), and [.github/run_tests.sh](../../../.github/run_tests.sh)
remain the authority.

Read with the [Induction gem hub](../induction.md), [HLL flux contract](hll-flux-contract.md),
[interpolation and staggering contract](interpolation-and-staggering-contract.md),
[gauge RHS contract](gauge-rhs-contract.md), and
[verification workflows](verification-workflows.md). Build, targeted run, and
CI routes belong in verification workflows, not here.

## Standalone Tests

Evidence labels used below distinguish direct replay source from configured
execution. `configure` discovers these test sources in both HDF5 modes;
`.github/run_tests.sh` and workflows list executions. A **compiled-unrun**
status requires a reviewer-observed binary; tracked source alone cannot prove
it.

| Behavior group | Test binary source | Required fixtures | Fixture owner |
| --- | --- | --- | --- |
| HLL flux with `B` and `Btilde` | [Unit_Tests/unit_test_HLL_flux.c](../../../Unit_Tests/unit_test_HLL_flux.c) | `HLL_flux_input.bin`, `HLL_flux_with_B_output.bin`, `HLL_flux_with_B_output_pert.bin`, `HLL_flux_with_Btilde_output.bin`, `HLL_flux_with_Btilde_output_pert.bin` | [Unit_Tests/data_gen/unit_test_data_HLL_flux.c](../../../Unit_Tests/data_gen/unit_test_data_HLL_flux.c) writes input and same-helper regression outputs. Speeds are random, not Flux_Source outputs. |
| ET Legacy HLL flux | [Unit_Tests/unit_test_ET_Legacy_HLL_flux.c](../../../Unit_Tests/unit_test_ET_Legacy_HLL_flux.c) | `ET_Legacy_HLL_flux_input.bin`, `ET_Legacy_HLL_flux_output.bin`, `ET_Legacy_HLL_flux_output_pert.bin` | [Unit_Tests/data_gen/unit_test_data_ET_Legacy_HLL_flux.c](../../../Unit_Tests/data_gen/unit_test_data_ET_Legacy_HLL_flux.c) writes input-side files; trusted outputs come through the external fixture download route visible in [.github/run_tests.sh](../../../.github/run_tests.sh). |
| Cell-centered ADM interpolation | [Unit_Tests/unit_test_induction_ccc_ADM.c](../../../Unit_Tests/unit_test_induction_ccc_ADM.c) | `induction_interpolation_input.bin`, `induction_interpolation_ADM_input.bin`, `induction_interpolation_ccc_ADM_output.bin`, `induction_interpolation_ccc_ADM_output_pert.bin` | [Unit_Tests/data_gen/unit_test_data_induction_interpolation.c](../../../Unit_Tests/data_gen/unit_test_data_induction_interpolation.c) writes shared input, ADM metric input, trusted output, and perturbation bars. |
| Cell-centered BSSN interpolation | [Unit_Tests/unit_test_induction_ccc_BSSN.c](../../../Unit_Tests/unit_test_induction_ccc_BSSN.c) | `induction_interpolation_input.bin`, `induction_interpolation_BSSN_input.bin`, `induction_interpolation_ccc_BSSN_output.bin`, `induction_interpolation_ccc_BSSN_output_pert.bin` | [Unit_Tests/data_gen/unit_test_data_induction_interpolation.c](../../../Unit_Tests/data_gen/unit_test_data_induction_interpolation.c) writes shared input, BSSN metric input, trusted output, and perturbation bars. |
| Vertex-centered ADM interpolation | [Unit_Tests/unit_test_induction_vvv_ADM.c](../../../Unit_Tests/unit_test_induction_vvv_ADM.c) | `induction_interpolation_input.bin`, `induction_interpolation_ADM_input.bin`, `induction_interpolation_vvv_ADM_output.bin`, `induction_interpolation_vvv_ADM_output_pert.bin` | [Unit_Tests/data_gen/unit_test_data_induction_interpolation.c](../../../Unit_Tests/data_gen/unit_test_data_induction_interpolation.c) writes shared input, ADM metric input, trusted output, and perturbation bars. |
| ET Legacy gauge RHS | [Unit_Tests/unit_test_ET_Legacy_induction_gauge_rhs.c](../../../Unit_Tests/unit_test_ET_Legacy_induction_gauge_rhs.c) | `ET_Legacy_induction_gauge_rhs_input.bin`, `ET_Legacy_induction_gauge_rhs_output.bin`, `ET_Legacy_induction_gauge_rhs_output_pert.bin` | [Unit_Tests/data_gen/unit_test_data_ET_Legacy_induction_gauge_rhs.c](../../../Unit_Tests/data_gen/unit_test_data_ET_Legacy_induction_gauge_rhs.c) writes input-side files; trusted outputs come through the external fixture download route visible in [.github/run_tests.sh](../../../.github/run_tests.sh). |

## Helper-Only Files

These `compute_*` files are helpers, not standalone tests:

- [Unit_Tests/compute_A_flux_with_B.c](../../../Unit_Tests/compute_A_flux_with_B.c)
  supports local and ET Legacy HLL flux replay using undensitized staggered `B`.
- [Unit_Tests/compute_A_flux_with_Btilde.c](../../../Unit_Tests/compute_A_flux_with_Btilde.c)
  supports local HLL flux replay using densitized `Btilde`.
- [Unit_Tests/compute_ccc_ADM.c](../../../Unit_Tests/compute_ccc_ADM.c)
  supports cell-centered ADM interpolation fixtures and replay.
- [Unit_Tests/compute_ccc_BSSN.c](../../../Unit_Tests/compute_ccc_BSSN.c)
  supports cell-centered BSSN interpolation fixtures and replay.
- [Unit_Tests/compute_vvv_ADM.c](../../../Unit_Tests/compute_vvv_ADM.c)
  supports vertex-centered ADM interpolation fixtures and replay.

## Data Generators

- [Unit_Tests/data_gen/unit_test_data_HLL_flux.c](../../../Unit_Tests/data_gen/unit_test_data_HLL_flux.c)
  writes `HLL_flux_input.bin`,
  `HLL_flux_with_B_output.bin`,
  `HLL_flux_with_B_output_pert.bin`,
  `HLL_flux_with_Btilde_output.bin`, and
  `HLL_flux_with_Btilde_output_pert.bin`.
- [Unit_Tests/data_gen/unit_test_data_induction_interpolation.c](../../../Unit_Tests/data_gen/unit_test_data_induction_interpolation.c)
  writes `induction_interpolation_input.bin`,
  `induction_interpolation_ADM_input.bin`,
  `induction_interpolation_BSSN_input.bin`,
  `induction_interpolation_ccc_ADM_output.bin`,
  `induction_interpolation_ccc_ADM_output_pert.bin`,
  `induction_interpolation_ccc_BSSN_output.bin`,
  `induction_interpolation_ccc_BSSN_output_pert.bin`,
  `induction_interpolation_vvv_ADM_output.bin`, and
  `induction_interpolation_vvv_ADM_output_pert.bin`.
- [Unit_Tests/data_gen/unit_test_data_ET_Legacy_HLL_flux.c](../../../Unit_Tests/data_gen/unit_test_data_ET_Legacy_HLL_flux.c)
  writes `ET_Legacy_HLL_flux_input.bin` and input perturbation data; the
  replay test reads `ET_Legacy_HLL_flux_output.bin` and
  `ET_Legacy_HLL_flux_output_pert.bin` from the downloaded ET Legacy fixtures.
- [Unit_Tests/data_gen/unit_test_data_ET_Legacy_induction_gauge_rhs.c](../../../Unit_Tests/data_gen/unit_test_data_ET_Legacy_induction_gauge_rhs.c)
  writes `ET_Legacy_induction_gauge_rhs_input.bin` and input perturbation data;
  the replay test reads `ET_Legacy_induction_gauge_rhs_output.bin` and
  `ET_Legacy_induction_gauge_rhs_output_pert.bin` from the downloaded ET Legacy
  fixtures.

## Behavior Notes

- HLL flux coverage compares `Ax_rhs`, `Ay_rhs`, and `Az_rhs` for both local
  `B` and `Btilde` variants. Contract detail belongs in
  [HLL flux contract](hll-flux-contract.md).
- HLL fixture generation chooses `cmin/cmax` independently from `[-10,10]`
  and rejects only exact zero sums. It does not validate production
  non-negative speed convention or Flux_Source coupling.
- ET Legacy HLL flux coverage replays legacy vector-potential HLL output for
  the `B` helper path.
- Interpolation coverage covers `ccc_ADM`, `ccc_BSSN`, and `vvv_ADM` fixture
  groups. Contract detail belongs in
  [interpolation and staggering contract](interpolation-and-staggering-contract.md).
- ET Legacy gauge RHS coverage assembles BSSN interpolation, vector-potential
  RHS terms, and `phitilde_rhs` replay. Contract detail belongs in
  [gauge RHS contract](gauge-rhs-contract.md).
- Induction tests listed here have no visible HDF5-only EOS table dependency.

Local HLL and interpolation reference outputs are produced through the same
public routines/helpers replayed by their tests. They detect regression against
stored results but are not independent implementations. ET Legacy output
fixtures are external comparisons; local ET data generators write input-side
files only.

## Configured Execution Status

- **Runner-listed:** `.github/run_tests.sh` invokes both ET Legacy tests, three
  interpolation tests, and local HLL test after downloading fixtures.
- **CI-configured:** all five compiler workflows contain `ET-Legacy`,
  `induction-interpolators`, and `induction-flux` jobs. Shared compile action
  builds tests and data generators; jobs execute downloaded fixtures.
- **Workflow-only for this change:** push/PR triggers ignore `wiki/**` and
  Markdown-only changes. Scheduled workflows still exist, but no execution was
  implied by tracked workflow configuration.

## Evidence Gaps

- No visible direct standalone `phitilde_rhs` test isolates only
  `ghl_calculate_phitilde_rhs`; the visible coverage is through
  [Unit_Tests/unit_test_ET_Legacy_induction_gauge_rhs.c](../../../Unit_Tests/unit_test_ET_Legacy_induction_gauge_rhs.c).
- No visible `vvv_BSSN` Induction test appears in the repo-local unit-test file
  list.
- `compute_*` helpers are not standalone tests.
- No HLL test feeds actual `ghl_calculate_characteristic_speed_dirn*` output
  into `ghl_HLL_vars`.
- No focused test covers HLL zero denominators, `Btilde` ET Legacy output, or
  vertex-wrapper `alpha`/`betai` remaining unwritten.

## Repo-Local References

- [wiki/test-map.md](../../test-map.md)
- [.github/run_tests.sh](../../../.github/run_tests.sh)
- [Unit_Tests/unit_test_HLL_flux.c](../../../Unit_Tests/unit_test_HLL_flux.c)
- [Unit_Tests/unit_test_ET_Legacy_HLL_flux.c](../../../Unit_Tests/unit_test_ET_Legacy_HLL_flux.c)
- [Unit_Tests/unit_test_induction_ccc_ADM.c](../../../Unit_Tests/unit_test_induction_ccc_ADM.c)
- [Unit_Tests/unit_test_induction_ccc_BSSN.c](../../../Unit_Tests/unit_test_induction_ccc_BSSN.c)
- [Unit_Tests/unit_test_induction_vvv_ADM.c](../../../Unit_Tests/unit_test_induction_vvv_ADM.c)
- [Unit_Tests/unit_test_ET_Legacy_induction_gauge_rhs.c](../../../Unit_Tests/unit_test_ET_Legacy_induction_gauge_rhs.c)
- [Unit_Tests/compute_A_flux_with_B.c](../../../Unit_Tests/compute_A_flux_with_B.c)
- [Unit_Tests/compute_A_flux_with_Btilde.c](../../../Unit_Tests/compute_A_flux_with_Btilde.c)
- [Unit_Tests/compute_ccc_ADM.c](../../../Unit_Tests/compute_ccc_ADM.c)
- [Unit_Tests/compute_ccc_BSSN.c](../../../Unit_Tests/compute_ccc_BSSN.c)
- [Unit_Tests/compute_vvv_ADM.c](../../../Unit_Tests/compute_vvv_ADM.c)
- [Unit_Tests/data_gen/unit_test_data_HLL_flux.c](../../../Unit_Tests/data_gen/unit_test_data_HLL_flux.c)
- [Unit_Tests/data_gen/unit_test_data_induction_interpolation.c](../../../Unit_Tests/data_gen/unit_test_data_induction_interpolation.c)
- [Unit_Tests/data_gen/unit_test_data_ET_Legacy_HLL_flux.c](../../../Unit_Tests/data_gen/unit_test_data_ET_Legacy_HLL_flux.c)
- [Unit_Tests/data_gen/unit_test_data_ET_Legacy_induction_gauge_rhs.c](../../../Unit_Tests/data_gen/unit_test_data_ET_Legacy_induction_gauge_rhs.c)
