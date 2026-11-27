# Reference Data And Test Drivers

> Perturbation-bar comparison rules (`ghl_pert_test_fail`), reference-data provenance classes, and binary data families across `Unit_Tests/`. · Status: contested
> Up: [Tests](index.md)

## Summary

Most `Unit_Tests/` replay drivers compare computed outputs against downloaded or generated trusted and perturbed reference binaries (`*_output.bin`/`*_output_pert.bin` or `*_unperturbed.bin`/`*_perturbed.bin`) via `ghl_pert_test_fail` and typed wrappers in `ghl_unit_tests.h` and `Unit_Tests/pert_test_fail_*.c`. Perturbation bars sample input-sensitivity tolerances rather than providing an independent physics oracle.

## Detail

### Comparison Tolerances And Nonfinite Rejection

- **Trusted, computed, and perturbed data:** "Trusted" denotes the stored comparison reference (downloaded from `GRHayL/TestData` or generated locally), not an independent analytic solution. Perturbed binaries establish per-point numerical sensitivity bars; exact binary `fread`/`fwrite` ordering lives in paired test and generator sources.
- **Strict nonfinite rejection:** `ghl_pert_test_fail` (`GRHayL/include/ghl_unit_tests.h`) rejects any nonfinite trusted, computed, or perturbed value.
- **Typed wrappers:** `pert_test_fail_conservatives.c` and `ghl_pert_test_fail_primitives_with_cutoffs` (`pert_test_fail_primitives.c`) apply small absolute roundoff allowances scaled to trusted momentum/velocity vectors (and Lorentz-factor conditioning for hybrid density recovery), whereas basic `ghl_pert_test_fail_primitives` compares componentwise without vector-scale allowances. `pert_test_fail_stress_energy.c` applies its absolute cutoff to the difference across the ten symmetric components.

Evidence (descriptive behavior): `GRHayL/include/ghl_unit_tests.h::ghl_pert_test_fail`; corroboration: `Unit_Tests/pert_test_fail_primitives.c::ghl_pert_test_fail_primitives_with_cutoffs`.

### Reference-Data Source Classes And Module Families

Claim status: contested; contradiction: CONTR-0003. See [CONTR-0003](../contradictions.md#contr-0003).

Claim status: contested; contradiction: CONTR-0004. See [CONTR-0004](../contradictions.md#contr-0004).

Reference data arrives via five distinct routes: local `Unit_Tests/data_gen/unit_test_data_*.c` binaries (compiled by `make datagen`, run manually), generation modes inside tests (`nrpyleakage_main.h` key `0`, `unit_test_con2prim_tabulated.c` key `0`), downloaded `GRHayL/TestData` `.bin` files, downloaded `stellarcollapse.org/EOS` `.h5.bz2` tables, and checked-in `Unit_Tests/sample_table/` files. Local PLM, WENO-z, Flux_Source, and Induction generators invoke the same library routines as replay, providing regression rather than independent-oracle checks.

| Module | Binary / table family | Blast-radius & routing notes |
| --- | --- | --- |
| Core | `grhayl_core_test_suite_input.bin` | Lowercase name shared across generator, test, runner, and workflows ([Core Tests](../core/tests-and-reference-data.md)). |
| Con2Prim / Core shared | `metric_Bfield_initial_data.bin`, `apply_conservative_limits_*`, `con2prim_multi_method_hybrid_*`, `enforce_primitive_limits_and_compute_u0_*`, `compute_conservs_and_Tmunu_*` | Editing `metric_Bfield_initial_data.bin` affects limits, recovery, `u0`, conservatives, and stress-energy tests. |
| Con2Prim tabulated | `con2prim_tabulated_*_{unperturbed,perturbed}.bin` | Paired with `SLy4` table and `pyghl append` ([Con2Prim Tests](../gems/con2prim/tests-and-reference-data.md)). |
| EOS | `simple_table.h5`, checked-in Hempel SFHo | Input tables, not perturbation bars ([HDF5 Sample Tables](hdf5-sample-tables.md)). |
| Flux_Source | `hybrid_flux_*`, `tabulated_flux_*` | Tabulated flux also reads `LS220` table ([Flux_Source Tests](../gems/flux-source/tests-and-reference-data.md)). |
| Induction | `induction_interpolation_*`, `HLL_flux_*` | Vector-potential HLL, distinct from Flux_Source HLLE ([Induction Tests](../gems/induction/tests-and-reference-data.md)). |
| Neutrinos | `nrpyleakage_{optically_thin_gas,constant_density_sphere,luminosities}_*` | Uses `*_unperturbed.bin`/`*_perturbed.bin` and `SLy4` ([Neutrinos Tests](../gems/neutrinos/tests-and-reference-data.md)). |
| Reconstruction & ET_Legacy | `PLM_reconstruction_*`, `WENOZ_reconstruction_*`, `ET_Legacy_*` | ET_Legacy outputs are externally supplied ([ET Legacy Comparison](et-legacy-comparison.md)). |

Evidence (descriptive behavior): `Unit_Tests/pert_test_fail_conservatives.c`; corroboration: `Unit_Tests/pert_test_fail_stress_energy.c`.

## Sources

- [GRHayL/include/ghl_unit_tests.h](../../GRHayL/include/ghl_unit_tests.h) - `ghl_pert_test_fail` and test helper declarations
- [Unit_Tests/pert_test_fail_conservatives.c](../../Unit_Tests/pert_test_fail_conservatives.c) - conservative comparator
- [Unit_Tests/pert_test_fail_primitives.c](../../Unit_Tests/pert_test_fail_primitives.c) - primitive comparator and `ghl_pert_test_fail_primitives_with_cutoffs`
- [Unit_Tests/pert_test_fail_stress_energy.c](../../Unit_Tests/pert_test_fail_stress_energy.c) - stress-energy comparator
- [.github/run_tests.sh](../../.github/run_tests.sh) - TestData and EOS table download steps

## See Also

- Depends on: [Data Generators And Helpers](data-generators-and-helpers.md)
- Depends on: [Runner And Build Outputs](runner-and-build-outputs.md)
- See also: [ET Legacy Comparison](et-legacy-comparison.md)
- See also: [Contradictions](../contradictions.md#contr-0004)
