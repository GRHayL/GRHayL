# Source Terms

> Geometric GRMHD source-term kernel (`ghl_calculate_source_terms`, `_checked`), caller-supplied metric derivatives, output fields, and EOS callback mutation. · Status: confirmed
> Up: [Flux_Source](index.md)

## Summary

`GRHayL/Flux_Source/ghl_calculate_source_terms.c` (generated from `generate_flux_source.py`) computes hydrodynamic spacetime-derivative source contributions for `cons->SD[0..2]` and `cons->tau` using caller-supplied metric spatial derivatives, extrinsic curvature, and global `ghl_compute_h_and_cs2`. It leaves `cons->rho`, `cons->Y_e`, and `cons->entropy` untouched and does not include neutrino matter sources.

## Detail

### Inputs, Output Fields, And EOS Callback Behavior

Declared in `GRHayL/include/ghl_flux_source.h` as `ghl_calculate_source_terms` (legacy `void`, aborts on EOS error) and `ghl_calculate_source_terms_checked` (returns `ghl_error_codes_t`), and built in both default and `--disable-hdf5` modes via `GRHayL/Flux_Source/make.code.defn`:
- **Caller inputs:** `eos`, mutable `prims` (with valid `u0`, rescaled `BU`, density, pressure, velocity), `metric` (ADM metric at the point), `metric_derivs_x`/`y`/`z` (caller-computed spatial derivatives of lapse, shift, and spatial metric packed into `ghl_metric_quantities`), `curv` (`ghl_extrinsic_curvature`), and `cons` output. Flux_Source never computes metric derivatives itself; callers maintain consistent centering and coordinate ordering.
- **Output & mutation:** Calls `ghl_compute_h_and_cs2(eos, prims, &h, &cs2)` (`cs2` unused). On callback failure, `_checked` returns the exact error and leaves `cons` untouched while any tabulated primitive clamping on `prims` is not rolled back (pass a copy if immutability is needed). On success, writes only `cons->SD[0..2]` and `cons->tau`; `rho`, `Y_e`, and `entropy` keep prior values.
- **Scope:** Computes stress-energy algebra inline without calling `ghl_compute_TDNmunu`/`TUPmunu`. Neutrino matter sources (`R_source`, `Q_source`) live separately in [Physics Model And Sources](../neutrinos/physics-model-and-sources.md).

Evidence (descriptive behavior): `ghl_calculate_source_terms.c::ghl_calculate_source_terms_checked`; corroboration: `ghl_flux_source.h::ghl_calculate_source_terms`.

### Test Coverage

- `Unit_Tests/unit_test_ET_Legacy_flux_source.c` (input generator `unit_test_data_ET_Legacy_flux_source.c`) exercises `ghl_calculate_source_terms` inside a combined flux-divergence/source replay using a test-local `compute_h_and_cs2` callback.
- `Unit_Tests/unit_test_tabulated_eos_compose.c` covers a production tabulated source call against analytic goldens in the Ubuntu-GCC `compose-regularized-eos` job.
- `Unit_Tests/unit_test_hybrid_flux.c` isolates checked source callback failure and verifies unchanged `cons` outputs.

## Sources

- [GRHayL/include/ghl_flux_source.h](../../../GRHayL/include/ghl_flux_source.h) - legacy and `_checked` source-term declarations
- [GRHayL/Flux_Source/ghl_calculate_source_terms.c](../../../GRHayL/Flux_Source/ghl_calculate_source_terms.c) - geometric source-term implementation
- [GRHayL/Flux_Source/generate_flux_source.py](../../../GRHayL/Flux_Source/generate_flux_source.py) - `generate_source()`
- [GRHayL/Flux_Source/make.code.defn](../../../GRHayL/Flux_Source/make.code.defn) - source-term build list
- [docs/raw/derivation.md](../../../docs/raw/derivation.md) - geometric source equations
- [docs/raw/Flux_Source.dox](../../../docs/raw/Flux_Source.dox) - Flux_Source documentation
- [Unit_Tests/unit_test_ET_Legacy_flux_source.c](../../../Unit_Tests/unit_test_ET_Legacy_flux_source.c) - ET_Legacy source replay
- [Unit_Tests/unit_test_tabulated_eos_compose.c](../../../Unit_Tests/unit_test_tabulated_eos_compose.c) - production tabulated source check

## See Also

- Depends on: [Variables And Conventions](../../physics/variables-and-conventions.md)
- Implements: [Generated NRPy Boundary](generated-nrpy-boundary.md)
- Validated by: [Flux_Source Tests And Reference Data](tests-and-reference-data.md)
- See also: [Physics Model And Sources](../neutrinos/physics-model-and-sources.md)
