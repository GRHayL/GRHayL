# Data Generators And Helpers

> `Unit_Tests/data_gen/unit_test_data_*.c` generators, shared C/H helper files, sanitizer notes, and sample-table pointers. · Status: confirmed
> Up: [Tests](index.md)

## Summary

`Unit_Tests/` pairs standalone `unit_test_*.c` drivers with `Unit_Tests/data_gen/unit_test_data_*.c` reference-data generators and helper-only C/H sources linked into tests rather than `libghl`. `make datagen` compiles generators under `test/data_gen/` without executing them, and `Unit_Tests/con2prim_newman_history.h` requires a bounds sanitizer to verify memory safety at final iterations.

## Detail

### Reference-Data Generators (`Unit_Tests/data_gen/`)

| Generator | Target subsystem | Generated data & limits |
| --- | --- | --- |
| `unit_test_data_grhayl_core_test_suite.c` | `GRHayL/GRHayL_Core/` | Writes lowercase `grhayl_core_test_suite_input.bin` via `ghl_randomize_metric`. (Ev: `Unit_Tests/data_gen/unit_test_data_grhayl_core_test_suite.c`) |
| `unit_test_data_con2prim_multi_method_hybrid.c` | `GRHayL/Con2Prim/Hybrid/` | Writes shared `metric_Bfield_initial_data.bin` plus hybrid C2P, limit, and stress-energy reference files. (Ev: `Unit_Tests/data_gen/unit_test_data_con2prim_multi_method_hybrid.c`) |
| `unit_test_data_hybrid_flux.c` / `unit_test_data_tabulated_flux.c` | `GRHayL/Flux_Source/` | Writes directional and entropy HLLE reference bars; tabulated requires `LS220_234r_136t_50y_analmu_20091212_SVNr26.h5` (HDF5-only). (Ev: `Unit_Tests/data_gen/unit_test_data_hybrid_flux.c`) |
| `unit_test_data_HLL_flux.c` / `unit_test_data_induction_interpolation.c` | `GRHayL/Induction/` | Writes vector-potential HLL (`B`/`Btilde`) and ADM/BSSN cell/vertex interpolation inputs and output bars. (Ev: `Unit_Tests/data_gen/unit_test_data_HLL_flux.c`) |
| `unit_test_data_PLM_reconstruction.c` / `unit_test_data_WENOZ_reconstruction.c` | `GRHayL/Reconstruction/` | Writes interior-stencil PLM (minmod/MC/superbee) and WENO-z reference bars with zero boundary placeholders; no PPM generator exists. (Ev: `Unit_Tests/data_gen/unit_test_data_PLM_reconstruction.c`) |
| `unit_test_data_ET_Legacy_{conservs,primitives,flux_source,HLL_flux,induction_gauge_rhs,reconstruction}.c` | ET_Legacy comparisons | Writes input and perturbed-input files with generator-specific random scales; trusted legacy outputs remain externally supplied. (Ev: `Unit_Tests/data_gen/unit_test_data_ET_Legacy_conservs.c`) |

### Shared Helper Files And Sanitizer Boundary

Helper-only files under `Unit_Tests/` and source-tree-only `GRHayL/include/ghl_unit_tests.h` are compiled into tests/generators, not installed into `libghl`:

| Helper file(s) | Role & defining test leaf |
| --- | --- |
| `con2prim_newman_history.h` | Drives tabulated Newman recovery through the final iteration; detecting an out-of-bounds history write requires running `unit_test_con2prim_tabulated` with a bounds sanitizer ([Con2Prim Tests](../gems/con2prim/tests-and-reference-data.md)). (Ev: `Unit_Tests/con2prim_newman_history.h`) |
| `ghl_test_helpers.c`, `ghl_test_helpers.h` | Shared Induction grid-header read, minimum-size, and overflow checks linked via `TOBJS` ([Induction Tests](../gems/induction/tests-and-reference-data.md)). (Ev: `Unit_Tests/ghl_test_helpers.h`) |
| `nrpyleakage_main.h` | Shared `main`/CLI key `0` (generate) and key `1` (replay) handler for table-backed NRPyLeakage tests; physics and classifier-fallback define `main(void)` ([Neutrinos Tests](../gems/neutrinos/tests-and-reference-data.md)). (Ev: `Unit_Tests/nrpyleakage_main.h`) |
| `pert_test_fail_{conservatives,primitives,stress_energy}.c` | Componentwise perturbation-bar comparators wrapping `ghl_pert_test_fail` ([Reference Data And Test Drivers](reference-data-and-test-drivers.md)). (Ev: `Unit_Tests/pert_test_fail_primitives.c`) |
| `randomize_metric.c`, `randomize_primitives.c` | Random metric (`lapse = 1`, `shift = 0`) and velocity/magnetic-field generators ([Core Tests](../core/tests-and-reference-data.md)). (Ev: `Unit_Tests/randomize_metric.c`) |
| `compute_A_flux_with_{B,Btilde}.c`, `compute_{ccc_ADM,ccc_BSSN,vvv_ADM}.c` | Induction HLL and ADM/BSSN interpolation test helpers ([Induction Tests](../gems/induction/tests-and-reference-data.md)). (Ev: `Unit_Tests/compute_A_flux_with_B.c`) |
| `tabulated_eos_unit_test_helpers.c`, `test_compute_h_and_cs2.c` | Analytic tabulated EOS oracle helpers vs. synthetic LS220-scale `ghl_compute_h_and_cs2` callback shared by tabulated flux tests (not a table oracle; [EOS Tests](../gems/eos/tests-and-reference-data.md)). (Ev: `Unit_Tests/test_compute_h_and_cs2.c`) |

For `Unit_Tests/sample_table/` (`Hempel_SFHoEOS_*_simple.h5`, `generate_simple_table.py`, `table_info.txt`), see [HDF5 Sample Tables](hdf5-sample-tables.md).

## Sources

- [Unit_Tests/data_gen/unit_test_data_grhayl_core_test_suite.c](../../Unit_Tests/data_gen/unit_test_data_grhayl_core_test_suite.c) - Core generator
- [Unit_Tests/data_gen/unit_test_data_con2prim_multi_method_hybrid.c](../../Unit_Tests/data_gen/unit_test_data_con2prim_multi_method_hybrid.c) - Con2Prim/Core shared generator
- [Unit_Tests/data_gen/unit_test_data_hybrid_flux.c](../../Unit_Tests/data_gen/unit_test_data_hybrid_flux.c) - hybrid flux generator
- [Unit_Tests/data_gen/unit_test_data_HLL_flux.c](../../Unit_Tests/data_gen/unit_test_data_HLL_flux.c) - Induction HLL generator
- [Unit_Tests/data_gen/unit_test_data_PLM_reconstruction.c](../../Unit_Tests/data_gen/unit_test_data_PLM_reconstruction.c) - PLM generator
- [Unit_Tests/data_gen/unit_test_data_ET_Legacy_conservs.c](../../Unit_Tests/data_gen/unit_test_data_ET_Legacy_conservs.c) - ET_Legacy conservs generator
- [Unit_Tests/con2prim_newman_history.h](../../Unit_Tests/con2prim_newman_history.h) - Newman history boundary test helper
- [Unit_Tests/ghl_test_helpers.h](../../Unit_Tests/ghl_test_helpers.h) - Induction grid-header helper declarations
- [Unit_Tests/nrpyleakage_main.h](../../Unit_Tests/nrpyleakage_main.h) - shared NRPyLeakage replay/generate driver
- [Unit_Tests/pert_test_fail_primitives.c](../../Unit_Tests/pert_test_fail_primitives.c) - primitive perturbation comparator
- [Unit_Tests/randomize_metric.c](../../Unit_Tests/randomize_metric.c) - metric randomizer
- [Unit_Tests/compute_A_flux_with_B.c](../../Unit_Tests/compute_A_flux_with_B.c) - Induction HLL helper
- [Unit_Tests/test_compute_h_and_cs2.c](../../Unit_Tests/test_compute_h_and_cs2.c) - synthetic flux EOS callback
- [GRHayL/include/ghl_unit_tests.h](../../GRHayL/include/ghl_unit_tests.h) - source-tree-only test helper header

## See Also

- Depends on: [Unit Test Catalog](unit-test-catalog.md)
- Depends on: [Reference Data And Test Drivers](reference-data-and-test-drivers.md)
- See also: [HDF5 Sample Tables](hdf5-sample-tables.md)
- See also: [ET Legacy Comparison](et-legacy-comparison.md)
