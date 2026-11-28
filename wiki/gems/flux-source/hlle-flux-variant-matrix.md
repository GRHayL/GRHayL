# HLLE Flux Variant Matrix

> Direct directional HLLE flux variants (`hybrid`, `hybrid_entropy`, `tabulated`, `tabulated_entropy`), wave-speed guards, EOS callback rules, and no-HDF5 behavior. · Status: confirmed
> Up: [Flux_Source](index.md)

## Summary

GRHayL exposes directional HLLE flux families (`hybrid`, `hybrid_entropy`, `tabulated`, `tabulated_entropy`) in `GRHayL/include/ghl_flux_source.h`, each with legacy `void` (abort-on-error) and `_checked` (`ghl_error_codes_t`) entry points across `dirn0`, `dirn1`, and `dirn2`. Deprecated unsuffixed `ghl_calculate_HLLE_fluxes_dirn0/1/2` function pointers in `ghl_eos_functions.h` are never assigned by GRHayL and have incompatible `const` signatures; new callers must call the direct suffixed functions.

## Detail

### Direct Variant Matrix And Output Fields

| Variant | Functions (`dirn0`, `dirn1`, `dirn2`; legacy & `_checked`) | Written `cons` fields | Manifest & replay test |
| --- | --- | --- | --- |
| `hybrid` | `ghl_calculate_HLLE_fluxes_dirn{0,1,2}_hybrid[_checked]` | `rho`, `tau`, `SD[0..2]` | `Flux_Source/hybrid/make.code.defn`; `unit_test_hybrid_flux.c`, `unit_test_ET_Legacy_flux_source.c` (Ev: `GRHayL/Flux_Source/hybrid/make.code.defn`) |
| `hybrid_entropy` | `ghl_calculate_HLLE_fluxes_dirn{0,1,2}_hybrid_entropy[_checked]` | `rho`, `tau`, `SD[0..2]`, `entropy` | `Flux_Source/hybrid_entropy/make.code.defn`; `unit_test_hybrid_flux.c` (Ev: `GRHayL/Flux_Source/hybrid_entropy/make.code.defn`) |
| `tabulated` | `ghl_calculate_HLLE_fluxes_dirn{0,1,2}_tabulated[_checked]` | `rho`, `tau`, `SD[0..2]`, `Y_e` | `Flux_Source/tabulated/make.code.defn`; `unit_test_tabulated_flux.c` (Ev: `GRHayL/Flux_Source/tabulated/make.code.defn`) |
| `tabulated_entropy` | `ghl_calculate_HLLE_fluxes_dirn{0,1,2}_tabulated_entropy[_checked]` | `rho`, `tau`, `SD[0..2]`, `Y_e`, `entropy` | `Flux_Source/tabulated_entropy/make.code.defn`; `unit_test_tabulated_flux.c` (Ev: `GRHayL/Flux_Source/tabulated_entropy/make.code.defn`) |

Only the listed fields of `cons` are written; unwritten fields retain prior values. Entropy variants read `prims_r->entropy` and `prims_l->entropy` directly.

### Wave-Speed Validation, EOS Callback, And `--disable-hdf5`

- **Wave-speed guards:** Negative `cmin`/`cmax` residues within `DBL_EPSILON * max(1, |cmin|, |cmax|)` are floored to `0`; larger negative values are rejected. Floored speeds must satisfy `cmin + cmax >= 1 / DBL_MAX` (finite) and non-overflowing product (within one rounding). Violations return `ghl_error_invalid_hlle_wavespeeds` from `_checked` before any EOS call or `cons` write (legacy wrappers abort). `_checked` does not validate geometry, `u0`, or output finiteness.
- **EOS callback & primitive mutation:** After wave-speed checks, kernels invoke global `ghl_compute_h_and_cs2` for `prims_r`, then `prims_l` (`cs2` is requested for single-callback compatibility with custom EOS integrations and unused here), returning the first callback error before writing `cons`. Production tabulated dispatch may clamp and mutate `prims_r`/`prims_l` in place without rollback on a second-call error; pass copies to preserve face states.
- **`--disable-hdf5` retention:** `configure --disable-hdf5` still compiles the real tabulated HLLE C kernels (excluding only tabulated flux tests/generators; Ubuntu-Clang's no-HDF5 `c2p-failure` job checks legacy tabulated symbol exports with `nm -D --defined-only` without linking a caller or inspecting `_checked` names). `README.md` states that these symbols are retained for link compatibility only, not for supported no-HDF5 execution: in that mode `ghl_initialize_eos_functions` cannot install table-backed `ghl_compute_h_and_cs2`, and the kernels call whatever the global pointer holds (null before EOS initialization, the hybrid routine after simple or hybrid initialization) without checking `eos->eos_type`.
- **Coverage limits:** Focused asymmetric analytic HLLE tests use zero magnetic field; no committed test verifies production tabulated callback mutation inside an HLLE call; deprecated generic HLLE pointers have no tests.

Evidence (descriptive behavior): `ghl_flux_source.h::ghl_calculate_HLLE_fluxes_dirn0_tabulated_checked`; corroboration: `generate_flux_source.py::generate_fluxes`.

Evidence (normative rule): `README.md::link compatibility`; corroboration: `initialize_eos.c::ghl_initialize_eos_functions`.

## Sources

- [GRHayL/include/ghl_flux_source.h](../../../GRHayL/include/ghl_flux_source.h) - direct HLLE legacy and `_checked` declarations
- [GRHayL/include/ghl_eos_functions.h](../../../GRHayL/include/ghl_eos_functions.h) - deprecated unassigned generic HLLE pointers
- [GRHayL/Flux_Source/hybrid/make.code.defn](../../../GRHayL/Flux_Source/hybrid/make.code.defn) - hybrid manifest
- [GRHayL/Flux_Source/hybrid_entropy/make.code.defn](../../../GRHayL/Flux_Source/hybrid_entropy/make.code.defn) - hybrid-entropy manifest
- [GRHayL/Flux_Source/tabulated/make.code.defn](../../../GRHayL/Flux_Source/tabulated/make.code.defn) - tabulated manifest
- [GRHayL/Flux_Source/tabulated_entropy/make.code.defn](../../../GRHayL/Flux_Source/tabulated_entropy/make.code.defn) - tabulated-entropy manifest
- [GRHayL/Flux_Source/generate_flux_source.py](../../../GRHayL/Flux_Source/generate_flux_source.py) - `generate_fluxes()`
- [GRHayL/GRHayL_Core/initialize_eos.c](../../../GRHayL/GRHayL_Core/initialize_eos.c) - no-HDF5 `ghl_compute_h_and_cs2` installation
- [README.md](../../../README.md) - tabulated HLLE link-compatibility rule for no-HDF5 builds
- [Unit_Tests/unit_test_hybrid_flux.c](../../../Unit_Tests/unit_test_hybrid_flux.c) - hybrid HLLE and wave-bound/error checks
- [Unit_Tests/unit_test_tabulated_flux.c](../../../Unit_Tests/unit_test_tabulated_flux.c) - tabulated HLLE checks

## See Also

- Depends on: [Characteristic Speeds](characteristic-speeds.md)
- Depends on: [Build And Configure](../../architecture/build-and-configure.md)
- Implements: [Generated NRPy Boundary](generated-nrpy-boundary.md)
- Validated by: [Flux_Source Tests And Reference Data](tests-and-reference-data.md)
