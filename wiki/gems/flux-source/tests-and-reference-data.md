# Flux_Source Tests And Reference Data

> Hybrid, tabulated, ET_Legacy, and CompOSE flux/source test drivers, data generators, synthetic callback limits, and Induction HLL exclusion. · Status: confirmed
> Up: [Flux_Source](index.md)

## Summary

Flux_Source validation spans `unit_test_hybrid_flux.c`, `unit_test_tabulated_flux.c` (LS220 + synthetic `test_compute_h_and_cs2.c`), `unit_test_ET_Legacy_flux_source.c`, and production-EOS checks in `unit_test_tabulated_eos_compose.c`. `unit_test_HLL_flux.c` and `unit_test_ET_Legacy_HLL_flux.c` test Induction vector-potential HLL fluxes, **not** Flux_Source hydrodynamic HLLE.

## Detail

### Direct Flux_Source Tests And Data Generators

| Test & generator | Data files | Scope & evidence boundaries |
| --- | --- | --- |
| `unit_test_hybrid_flux.c` / `data_gen/unit_test_data_hybrid_flux.c` | `hybrid_flux_{input,output,output_pert}.bin` | Replays `hybrid` and `hybrid_entropy` (all directions) using stored wave speeds, plus wave-bound/error checks, checked speed/source callback fault injection, and independent zero-`B` asymmetric HLLE algebra. Generator perturbs up to `1e-14` and rescales velocities to `W in [1.2, 5]` (below `W_max = 10`). (Ev: `unit_test_hybrid_flux.c`) |
| `unit_test_tabulated_flux.c` / `data_gen/unit_test_data_tabulated_flux.c` | `tabulated_flux_{input,output,output_pert}.bin`, `LS220_234r_136t_50y_analmu_20091212_SVNr26.h5` | HDF5-only; replays `tabulated` and `tabulated_entropy` (all directions) plus wave-bound/error and zero-`B` checks. Replaces production EOS dispatch with synthetic `ghl_test_compute_h_and_cs2` (never reads LS220 sound speed). Generator perturbs up to `1e-12` with `W in [1.2, 5]`. (Ev: `unit_test_tabulated_flux.c`) |
| `unit_test_ET_Legacy_flux_source.c` / `data_gen/unit_test_data_ET_Legacy_flux_source.c` | `ET_Legacy_flux_source_{input,output,output_pert}.bin` | Recomputes characteristic speeds (with a test-local `compute_h_and_cs2` callback), hybrid HLLE fluxes, flux divergence, and source terms on equal spacings `dx = dy = dz = 0.1`; outputs are external legacy data. (Ev: `unit_test_ET_Legacy_flux_source.c`) |
| `unit_test_tabulated_eos_compose.c` | Converted CompOSE HDF5 table (CI Ubuntu-GCC `compose-regularized-eos` only) | Checks analytic magnetized characteristic speeds, production tabulated table-bound/mutation behavior, identical-state HLLE/entropy fluxes, and production source terms. (Ev: `unit_test_tabulated_eos_compose.c`) |

### Induction Exclusion, Runner Selection, And Coverage Gaps

- **Not Flux_Source coverage:** `unit_test_HLL_flux.c`, `unit_test_ET_Legacy_HLL_flux.c`, and `unit_test_data_HLL_flux.c` belong to [Induction Tests And Reference Data](../induction/tests-and-reference-data.md).
- **Runner & CI:** `.github/run_tests.sh` downloads Flux data via `.github/et-legacy-testdata-ref` plus `LS220` and runs ET_Legacy, hybrid, and tabulated flux tests; compiler workflows run `flux` (`hybrid_flux`, `tabulated_flux`—downloading LS220 only on the tabulated leg) and `ET-Legacy` (`flux_source`). `--disable-hdf5` excludes `unit_test_tabulated_flux` and its generator.
- **Coverage gaps:** No independent analytic oracle for magnetized HLLE terms; ET_Legacy uses isotropic spacings `dx = dy = dz = 0.1`; production tabulated callback mutation is checked for speeds, not HLLE; no test injects failure into the production tabulated speed callback; no test verifies unwritten `cons` fields stay untouched on success; deprecated generic HLLE pointers are untested; no test connects Flux_Source speed output to Induction HLL input.

Evidence (descriptive behavior): `test_compute_h_and_cs2.c`; corroboration: `run_tests.sh::unit_test_tabulated_flux`.

## Sources

- [Unit_Tests/unit_test_hybrid_flux.c](../../../Unit_Tests/unit_test_hybrid_flux.c) - hybrid HLLE and callback-fault tests
- [Unit_Tests/data_gen/unit_test_data_hybrid_flux.c](../../../Unit_Tests/data_gen/unit_test_data_hybrid_flux.c) - hybrid flux generator
- [Unit_Tests/unit_test_tabulated_flux.c](../../../Unit_Tests/unit_test_tabulated_flux.c) - tabulated HLLE test
- [Unit_Tests/data_gen/unit_test_data_tabulated_flux.c](../../../Unit_Tests/data_gen/unit_test_data_tabulated_flux.c) - tabulated flux generator
- [Unit_Tests/unit_test_ET_Legacy_flux_source.c](../../../Unit_Tests/unit_test_ET_Legacy_flux_source.c) - ET_Legacy flux/source test
- [Unit_Tests/data_gen/unit_test_data_ET_Legacy_flux_source.c](../../../Unit_Tests/data_gen/unit_test_data_ET_Legacy_flux_source.c) - ET_Legacy input generator
- [Unit_Tests/unit_test_tabulated_eos_compose.c](../../../Unit_Tests/unit_test_tabulated_eos_compose.c) - production tabulated speed/flux/source checks
- [Unit_Tests/test_compute_h_and_cs2.c](../../../Unit_Tests/test_compute_h_and_cs2.c) - synthetic LS220-scale callback
- [.github/run_tests.sh](../../../.github/run_tests.sh) - flux downloads and execution order

## See Also

- Depends on: [HLLE Flux Variant Matrix](hlle-flux-variant-matrix.md)
- Depends on: [Characteristic Speeds](characteristic-speeds.md)
- Contrasts with: [Induction Tests And Reference Data](../induction/tests-and-reference-data.md)
- See also: [Reference Data And Test Drivers](../../tests/reference-data-and-test-drivers.md)
