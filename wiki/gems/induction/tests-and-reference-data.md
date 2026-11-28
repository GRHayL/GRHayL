# Induction Tests And Reference Data

> Induction unit-test binaries, private packing/header helpers, reference-data generators, stencil size guards, and documented coverage gaps. · Status: confirmed
> Up: [Induction](index.md)

## Summary

Induction coverage consists of HDF5-independent C tests (`unit_test_HLL_flux.c`, `unit_test_ET_Legacy_HLL_flux.c`, `unit_test_induction_{ccc_ADM,ccc_BSSN,vvv_ADM}.c`, `unit_test_ET_Legacy_induction_gauge_rhs.c`) supported by private `compute_*` packing helpers and `ghl_test_helpers.{c,h}` linked via `TOBJS`. Local HLL and interpolation generators call the same helpers as replay; no standalone `phitilde_rhs` test, `vvv_BSSN` test, or Flux_Source-to-Induction coupling test exists.

## Detail

### Test Binaries, Helpers, And Reference-Data Generators

| Test binary | Reference files | Generator & helper notes |
| --- | --- | --- |
| `unit_test_HLL_flux.c` | `HLL_flux_input.bin`, `HLL_flux_with_{B,Btilde}_output{,_pert}.bin` | `unit_test_data_HLL_flux.c` + `compute_A_flux_with_{B,Btilde}.c`; random speeds in `[-10, 10]` rejecting only exact zero sums; zero-initialized `index = 0` planes; requires grid `>= 5`. (Ev: `unit_test_HLL_flux.c`) |
| `unit_test_ET_Legacy_HLL_flux.c` | `ET_Legacy_HLL_flux_{input,output,output_pert}.bin` | `unit_test_data_ET_Legacy_HLL_flux.c` (inputs only) + `compute_A_flux_with_B.c`; external legacy outputs cover `B` path only. (Ev: `unit_test_ET_Legacy_HLL_flux.c`) |
| `unit_test_induction_{ccc_ADM,ccc_BSSN,vvv_ADM}.c` | `induction_interpolation_{input,ADM_input,BSSN_input,{ccc_ADM,ccc_BSSN,vvv_ADM}_output{,_pert}}.bin` | `unit_test_data_induction_interpolation.c` + `compute_{ccc_ADM,ccc_BSSN,vvv_ADM}.c`; requires grid `>= 3`; `vvv_ADM` includes an independent asymmetric affine test checking orientation and untouched `alpha`/`betai`. (Ev: `unit_test_induction_vvv_ADM.c`) |
| `unit_test_ET_Legacy_induction_gauge_rhs.c` | `ET_Legacy_induction_gauge_rhs_{input,output,output_pert}.bin` | `unit_test_data_ET_Legacy_induction_gauge_rhs.c` (inputs only) + `compute_ccc_BSSN.c`; requires grid `>= 7`. (Ev: `unit_test_ET_Legacy_induction_gauge_rhs.c`) |

`Unit_Tests/ghl_test_helpers.c` and `.h` (private helpers linked through `TOBJS`) enforce grid-header reads, minimum stencil lengths (`5` for HLL, `3` for interpolation, `7` for ET_Legacy gauge RHS), and overflow-safe cube allocations across every Induction replay binary.

### Coverage Gaps

- No standalone `ghl_calculate_phitilde_rhs`-only test (only combined BSSN ET_Legacy replay).
- No `vvv_BSSN` declaration, helper, or test.
- No test feeding actual `ghl_calculate_characteristic_speed_dirn*` outputs into `ghl_HLL_vars`, no zero-denominator HLL test, and no `Btilde` ET_Legacy comparison.

Evidence (descriptive behavior): `ghl_test_helpers.c`; corroboration: `unit_test_data_induction_interpolation.c`.

## Sources

- [Unit_Tests/unit_test_HLL_flux.c](../../../Unit_Tests/unit_test_HLL_flux.c) - local `B` and `Btilde` HLL replay
- [Unit_Tests/unit_test_ET_Legacy_HLL_flux.c](../../../Unit_Tests/unit_test_ET_Legacy_HLL_flux.c) - ET_Legacy `B` HLL replay
- [Unit_Tests/unit_test_induction_ccc_ADM.c](../../../Unit_Tests/unit_test_induction_ccc_ADM.c) - cell-centered ADM test
- [Unit_Tests/unit_test_induction_ccc_BSSN.c](../../../Unit_Tests/unit_test_induction_ccc_BSSN.c) - cell-centered BSSN test
- [Unit_Tests/unit_test_induction_vvv_ADM.c](../../../Unit_Tests/unit_test_induction_vvv_ADM.c) - vertex-centered ADM test + affine check
- [Unit_Tests/unit_test_ET_Legacy_induction_gauge_rhs.c](../../../Unit_Tests/unit_test_ET_Legacy_induction_gauge_rhs.c) - ET_Legacy gauge RHS test
- [Unit_Tests/ghl_test_helpers.c](../../../Unit_Tests/ghl_test_helpers.c) - shared grid-header and minimum-size validator
- [Unit_Tests/data_gen/unit_test_data_induction_interpolation.c](../../../Unit_Tests/data_gen/unit_test_data_induction_interpolation.c) - interpolation reference generator
- [Unit_Tests/data_gen/unit_test_data_HLL_flux.c](../../../Unit_Tests/data_gen/unit_test_data_HLL_flux.c) - HLL reference generator

## See Also

- Depends on: [HLL Flux](hll-flux.md)
- Depends on: [Interpolation And Staggering](interpolation-and-staggering.md)
- Depends on: [Gauge RHS](gauge-rhs.md)
- See also: [Induction Verification Workflows](verification-workflows.md)
