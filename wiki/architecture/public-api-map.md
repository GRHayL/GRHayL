# Public API Map

> Module mapping, installation status, header aggregation, and downstream impact for every header in `GRHayL/include/`. · Status: confirmed
> Up: [Architecture](index.md)

## Summary

Every public header under `GRHayL/include/` declares candidate functions, structs, or enums; installed availability comes from `GRHayL/include/make.code.defn` (`install_headers`), linkability requires a built definition or storage, and function-pointer callability requires compatible runtime initialization. `ghl_unit_tests.h` and `ghl_eos_functions_declaration.h` are not installed.

## Detail

### Core, Debug, IO, And Metric Headers

| Header | Module & install status | Implementations, docs & tests | Downstream impact |
| --- | --- | --- | --- |
| [ghl.h](../../GRHayL/include/ghl.h) | Core shared header (installed); includes `ghl_io.h`, `ghl_metric_helpers.h`, EOS function pointers, and `ghl_debug.h`. (Ev: `GRHayL/include/ghl.h`) | `GRHayL/GRHayL_Core/`, `docs/raw/GRHayL_Core.dox`, `unit_test_grhayl_core_test_suite.c`, `unit_test_compute_conservs_and_Tmunu.c`. | Highest risk: `GRHayLib` `initialize_and_shutdown.c`, `param.ccl`, and module headers depend on its structs/enums. |
| [ghl_io.h](../../GRHayL/include/ghl_io.h) | Core logging/error API (installed). (Ev: `GRHayL/include/ghl_io.h`) | `info_warn_error.c`, `abort_if_error.c`, `docs/raw/GRHayL_Core.dox`, `unit_test_code_error.c`. | Macro or exit-code changes affect every module and downstream Cactus error handling. |
| [ghl_metric_helpers.h](../../GRHayL/include/ghl_metric_helpers.h) | Core header-only metric vector helpers (installed). (Ev: `GRHayL/include/ghl_metric_helpers.h`) | Header-only; paired with `compute_ADM_auxiliaries.c`, `unit_test_grhayl_core_test_suite.c`. | Indexing layout feeds stress-energy, Con2Prim, Flux_Source, and Induction. |
| [ghl_debug.h](../../GRHayL/include/ghl_debug.h) | Shared header-only primitive/conservative print helpers (installed). (Ev: `GRHayL/include/ghl_debug.h`) | Uses `info_warn_error.c`; exercised (without output assertions) by `unit_test_con2prim_debug.c`. | Formatting edits affect every caller-compiled debug dump. |
| [ghl_eos_functions.h](../../GRHayL/include/ghl_eos_functions.h) | Shared EOS extern function-pointer declarations (installed). (Ev: `GRHayL/include/ghl_eos_functions.h`) | `GRHayL/EOS/{Hybrid,Tabulated}/`, `initialize_eos.c`, `docs/raw/EOS.dox`, EOS/flux/C2P tests. | Global dispatch must match active `ghl_eos_parameters`; deprecated `ghl_calculate_HLLE_fluxes_dirn*` pointers stay unassigned. |
| [ghl_eos_functions_declaration.h](../../GRHayL/include/ghl_eos_functions_declaration.h) | Internal storage definitions for EOS/HLLE function pointers (**not installed**). (Ev: `GRHayL/include/ghl_eos_functions_declaration.h`) | Included by `initialize_eos.c` and EOS implementations; omitted from `install_headers`. | Declaration/definition mismatch breaks configured builds and direct-compile thorns. |
| [ghl_unit_tests.h](../../GRHayL/include/ghl_unit_tests.h) | Source-tree test-helper aggregate (**not installed**, excluded by `Doxyfile`). (Ev: `GRHayL/include/ghl_unit_tests.h`) | Inline helpers plus `Unit_Tests/` and `Unit_Tests/data_gen/` definitions. | Internal test API only; external callers must not treat it as installed `libghl` API. |

For the `GRHayL_Core/make.code.defn` source list and the non-Core `ghl.h` declaration audit (`ghl_get_con2prim_routine_name`, `ghl_compute_SU_Bsq_Ssq_BdotS`), see [Shared Parameters And Enums](../core/shared-parameters-and-enums.md). Struct comments in `ghl.h` reference an external GitHub wiki URL (`https://github.com/GRHayL/GRHayL/wiki/<struct>`); repo-local headers, `docs/raw/`, and `wiki/` are authoritative.

Evidence (descriptive behavior): `GRHayL/include/make.code.defn::install_headers`; corroboration: `GRHayL/include/ghl.h`.

### Gem Public Headers

| Header | Module & install status | Implementations, docs & tests | Key availability & downstream notes |
| --- | --- | --- | --- |
| [ghl_atmosphere.h](../../GRHayL/include/ghl_atmosphere.h) | Atmosphere (installed). (Ev: `GRHayL/include/ghl_atmosphere.h`) | `set_prims_to_constant_atm.c`, `docs/raw/Atmosphere.dox`, `unit_test_grhayl_core_test_suite.c`, `unit_test_c2p_nn_guess.c`. | Included by `GRHayLib.h`; `ghl_set_prims_to_radial_falloff_atm` is declared but unbuilt. |
| [ghl_con2prim.h](../../GRHayL/include/ghl_con2prim.h) | Con2Prim (installed). (Ev: `GRHayL/include/ghl_con2prim.h`) | `GRHayL/Con2Prim/{,Hybrid,Tabulated}/`, `docs/raw/Con2Prim.dox`, C2P unit tests. | `Noble1D_entropy2` requires one-piece constant-Gamma; without HDF5, direct tabulated solvers become non-mutating stubs. |
| [ghl_nrpyeos_hybrid.h](../../GRHayL/include/ghl_nrpyeos_hybrid.h) | EOS hybrid/piecewise polytrope (installed). (Ev: `GRHayL/include/ghl_nrpyeos_hybrid.h`) | `GRHayL/EOS/Hybrid/`, `docs/raw/EOS.dox`, `unit_test_piecewise_polytrope.c`. | Included by `GRHayLib.h`; affects simple/hybrid initialization and Con2Prim compatibility. |
| [ghl_nrpyeos_tabulated.h](../../GRHayL/include/ghl_nrpyeos_tabulated.h) | EOS tabulated HDF5 API (installed). (Ev: `GRHayL/include/ghl_nrpyeos_tabulated.h`) | `GRHayL/EOS/Tabulated/`, `docs/raw/EOS.dox`, `unit_test_tabulated_eos.c`, C2P/Neutrinos tests. | Must agree with `GHL_DISABLE_HDF5`, table memory cleanup, and `GRHayLib` `requires HDF5`. |
| [ghl_flux_source.h](../../GRHayL/include/ghl_flux_source.h) | Flux_Source (installed). (Ev: `GRHayL/include/ghl_flux_source.h`) | `GRHayL/Flux_Source/` + 4 HLLE subdirs, `docs/raw/Flux_Source.dox`, hybrid/tabulated/ET_Legacy flux tests. | Legacy `void` names abort on error; `_checked` variants return `ghl_error_codes_t`; direct tabulated HLLE stays built without HDF5. |
| [ghl_induction.h](../../GRHayL/include/ghl_induction.h) | Induction (installed). (Ev: `GRHayL/include/ghl_induction.h`) | `GRHayL/Induction/` + `Interpolators/`, `docs/raw/Induction.dox`, HLL/ccc/vvv/gauge tests. | Stencil and staggering conventions are caller-managed; `GRHayLib.h` aggregates it. |
| [ghl_radiation.h](../../GRHayL/include/ghl_radiation.h) & [ghl_nrpyleakage.h](../../GRHayL/include/ghl_nrpyleakage.h) | Neutrinos / NRPyLeakage (both installed). (Ev: `GRHayL/include/ghl_radiation.h`) | `GRHayL/Neutrinos/NRPyLeakage/`, inline Doxygen group in `ghl_radiation.h` (no `.dox` file), NRPyLeakage tests. | `ghl_radiation.h` defines structs then includes `ghl_nrpyleakage.h`; table-backed tests couple to EOS/HDF5. |
| [ghl_reconstruction.h](../../GRHayL/include/ghl_reconstruction.h) | Reconstruction (installed). (Ev: `GRHayL/include/ghl_reconstruction.h`) | `GRHayL/Reconstruction/{PLM,PPM,WENOZ}/`, `docs/raw/Reconstruction.dox` (read-only), PLM/WENOZ/ET_Legacy tests. | Built public API with no in-repo production caller; caller supplies stencils and face conventions. |

Claim status: contested; contradiction: CONTR-0009. See [CONTR-0009](../contradictions.md#contr-0009).

### Header Aggregation And GRHayLib Boundary

`implementations/GRHayLib/src/GRHayLib.h` aggregates GRHayL's installed module headers and declares the thorn globals `ghl_params` and `ghl_eos`. Any public header rename, split, struct/enum change (`ghl_con2prim_id_t`, `ghl_eos_parameters`, `ghl_parameters`), or EOS/Con2Prim initializer signature change requires synchronized review of `implementations/GRHayLib/param.ccl`, `src/GRHayLib.h`, and `src/initialize_and_shutdown.c` per [GRHayLib Runtime Parameters](../implementations/grhaylib-runtime-parameters.md).

Evidence (descriptive behavior): `implementations/GRHayLib/src/GRHayLib.h::ghl_params`; corroboration: `GRHayL/include/ghl_radiation.h::ghl_nrpyleakage.h`.

## Sources

- [GRHayL/include/make.code.defn](../../GRHayL/include/make.code.defn) - `install_headers` list
- [GRHayL/include/ghl.h](../../GRHayL/include/ghl.h) - Core structs, enums, and transitive includes
- [GRHayL/include/ghl_io.h](../../GRHayL/include/ghl_io.h) - error and logging declarations
- [GRHayL/include/ghl_metric_helpers.h](../../GRHayL/include/ghl_metric_helpers.h) - inline metric vector helpers
- [GRHayL/include/ghl_debug.h](../../GRHayL/include/ghl_debug.h) - inline debug print helpers
- [GRHayL/include/ghl_eos_functions.h](../../GRHayL/include/ghl_eos_functions.h) - EOS and deprecated HLLE extern pointers
- [GRHayL/include/ghl_eos_functions_declaration.h](../../GRHayL/include/ghl_eos_functions_declaration.h) - internal pointer storage definitions
- [GRHayL/include/ghl_unit_tests.h](../../GRHayL/include/ghl_unit_tests.h) - source-tree test helper aggregate
- [GRHayL/include/ghl_atmosphere.h](../../GRHayL/include/ghl_atmosphere.h) - Atmosphere declarations
- [GRHayL/include/ghl_con2prim.h](../../GRHayL/include/ghl_con2prim.h) - Con2Prim declarations
- [GRHayL/include/ghl_nrpyeos_hybrid.h](../../GRHayL/include/ghl_nrpyeos_hybrid.h) - Hybrid EOS declarations
- [GRHayL/include/ghl_nrpyeos_tabulated.h](../../GRHayL/include/ghl_nrpyeos_tabulated.h) - Tabulated EOS declarations
- [GRHayL/include/ghl_flux_source.h](../../GRHayL/include/ghl_flux_source.h) - Flux_Source declarations
- [GRHayL/include/ghl_induction.h](../../GRHayL/include/ghl_induction.h) - Induction declarations
- [GRHayL/include/ghl_radiation.h](../../GRHayL/include/ghl_radiation.h) - radiation structs and Neutrinos Doxygen group
- [GRHayL/include/ghl_nrpyleakage.h](../../GRHayL/include/ghl_nrpyleakage.h) - NRPyLeakage declarations
- [GRHayL/include/ghl_reconstruction.h](../../GRHayL/include/ghl_reconstruction.h) - Reconstruction declarations
- [implementations/GRHayLib/src/GRHayLib.h](../../implementations/GRHayLib/src/GRHayLib.h) - GRHayLib aggregate header

## See Also

- Depends on: [Shared Parameters And Enums](../core/shared-parameters-and-enums.md)
- Depends on: [Build And Configure](build-and-configure.md)
- See also: [GRHayLib Runtime Parameters](../implementations/grhaylib-runtime-parameters.md)
- See also: [Change Impact](change-impact.md)
