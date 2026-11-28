# GRHayLib Cactus Build

> Thorn overview, Cactus CCL configuration, aggregate header, direct-compile symlink manifest, and host-only CarpetX caveat for `implementations/GRHayLib/`. · Status: confirmed
> Up: [Implementations](index.md)

## Summary

`implementations/GRHayLib/` is a Cactus/Einstein Toolkit thorn that directly compiles GRHayL source directories inside a Cactus thorn build rather than linking an externally installed `libghl`. It declares `PROVIDES GRHayL`, hard-codes `requires HDF5`, exports `GRHayLib.h` (aggregating installed GRHayL headers plus global pointers `ghl_params` and `ghl_eos`), and currently provides host-only functions for CarpetX simulations.

## Detail

### Thorn Purpose And Cactus CCL Declarations

`implementations/GRHayLib/README` describes GRHayLib as providing GRHayL core features and initializing the shared parameter and EOS structs used across downstream thorns, while noting that CarpetX support is currently host-only (device execution remains unfinished). Static CCL declarations establish the thorn build interface:

- **`configuration.ccl`:** declares `PROVIDES GRHayL` and hard-codes `requires HDF5`. Unlike standalone GRHayL `./configure --disable-hdf5`, the Cactus thorn provides no no-HDF5 build mode.
- **`interface.ccl`:** declares `implements: GRHayLib` and `INCLUDE HEADER: GRHayLib.h in GRHayLib.h`.
- **`schedule.ccl`:** schedules `GRHayLib_initialize` at `CCTK_WRAGH` only when `ID_converter_ILGRMHD` is not active, and always schedules `GRHayLib_terminate` at `CCTK_TERMINATE`.

Evidence (descriptive behavior): `configuration.ccl::requires HDF5`; corroboration: `interface.ccl::GRHayLib.h` and `schedule.ccl::CCTK_WRAGH`.

### Aggregate Header And Direct-Compile Source Symlinks

- **`src/GRHayLib.h`:** aggregates GRHayL's installed public module headers and declares `extern ghl_eos_parameters *ghl_eos` and `extern ghl_parameters *ghl_params`. Static header parity requires every header directly included by `GRHayLib.h` to exist in `GRHayL/include/` and appear in `GRHayL/include/make.code.defn` (`install_headers`); the tracked symlink `implementations/GRHayLib/src/include` points directly to `GRHayL/include/`.
- **`src/make.code.defn`:** lists thorn-local `SRCS = initialize_and_shutdown.c` and `SUBDIRS` covering every upstream GRHayL source-bearing module directory (including `Con2Prim/Tabulated/neural_network_guess`, plus intermediate `Con2Prim/Hybrid`), exposed via tracked directory symlinks under `implementations/GRHayLib/src/`.

Static CCL, header, and symlink parity can be checked in this repository, whereas Cactus compilation, scheduling, parameter parsing, and runtime execution require an external Cactus/Einstein Toolkit environment; upstream `Unit_Tests/` and `ET_Legacy` tests do not verify the thorn.

Evidence (descriptive behavior): `implementations/GRHayLib/src/make.code.defn::initialize_and_shutdown.c`; corroboration: `GRHayLib.h::ghl_eos` and `README::CarpetX`.

## Sources

- [implementations/GRHayLib/README](../../implementations/GRHayLib/README) - thorn overview and host-only CarpetX note
- [implementations/GRHayLib/doc/documentation.tex](../../implementations/GRHayLib/doc/documentation.tex) - ThornGuide documentation
- [implementations/GRHayLib/configuration.ccl](../../implementations/GRHayLib/configuration.ccl) - `PROVIDES GRHayL`, `requires HDF5`
- [implementations/GRHayLib/interface.ccl](../../implementations/GRHayLib/interface.ccl) - `implements: GRHayLib`, exported `GRHayLib.h`
- [implementations/GRHayLib/schedule.ccl](../../implementations/GRHayLib/schedule.ccl) - `CCTK_WRAGH` and `CCTK_TERMINATE` scheduling
- [implementations/GRHayLib/src/GRHayLib.h](../../implementations/GRHayLib/src/GRHayLib.h) - aggregate header and `ghl_params`/`ghl_eos` declarations
- [implementations/GRHayLib/src/make.code.defn](../../implementations/GRHayLib/src/make.code.defn) - `SRCS` and direct-compile `SUBDIRS`
- [GRHayL/include/make.code.defn](../../GRHayL/include/make.code.defn) - upstream `install_headers` parity reference

## See Also

- Depends on: [Public API Map](../architecture/public-api-map.md)
- Contrasts with: [Build And Configure](../architecture/build-and-configure.md)
- See also: [GRHayLib Runtime Parameters](grhaylib-runtime-parameters.md)
- See also: [GRHayLib Verification](grhaylib-verification.md)
