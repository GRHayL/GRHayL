# GRHayLib Hub

This page routes GRHayLib KB work for the Cactus/Einstein Toolkit thorn that
directly compiles GRHayL code. It defines the ownership boundary only; route
physics and gem behavior to existing Core and gem pages.

Repo-local files under [`implementations/GRHayLib/`](../../implementations/GRHayLib/)
remain the ground truth for this thorn.

## Ownership Boundary

GRHayLib is a Cactus/Einstein Toolkit thorn for direct GRHayL compilation and
direct-compile access to GRHayL. It differs from normal installed-library use
by letting the downstream infrastructure compile GRHayL code inside the thorn
build.

Public thorn access is through [`src/GRHayLib.h`](../../implementations/GRHayLib/src/GRHayLib.h).
That header aggregates GRHayL module headers and exposes the global pointers
`ghl_params` and `ghl_eos`. The globals are owned by GRHayLib lifecycle code in
[`src/initialize_and_shutdown.c`](../../implementations/GRHayLib/src/initialize_and_shutdown.c):
`GRHayLib_initialize` allocates and initializes them from Cactus parameters, and
`GRHayLib_terminate` frees them, including `ghl_tabulated_free_memory` for
tabulated EOS state.

That ownership statement describes the normal initialized path. Allocation is
unchecked, the schedule can skip initialization when `ID_converter_ILGRMHD` is
active while still scheduling termination, and no alternate owner is present
in this checkout. Use the runtime contract before assuming globals are valid.

The implementation README says GRHayLib provides core library features and
initializes the GRHayL parameter and EOS structs used through the simulation.
It also records a CarpetX caveat: GRHayLib currently provides functions on the
host only for CarpetX simulations; host and device access remains ongoing work.

## Ground Truth Files

- [`README`](../../implementations/GRHayLib/README): thorn purpose,
  downstream thorns, and host-only CarpetX caveat.
- [`doc/documentation.tex`](../../implementations/GRHayLib/doc/documentation.tex):
  ThornGuide source for `GRHayLib.h`, global pointers, EOS setup, Con2Prim
  backups, primitive guesses, and steering caveats.
- [`configuration.ccl`](../../implementations/GRHayLib/configuration.ccl):
  `PROVIDES GRHayL` and hard `requires HDF5`.
- [`interface.ccl`](../../implementations/GRHayLib/interface.ccl):
  `implements: GRHayLib` and exported `GRHayLib.h`.
- [`param.ccl`](../../implementations/GRHayLib/param.ccl): Cactus-facing
  GRHayL, EOS, Con2Prim, limits, entropy/temperature, and PPM parameters.
- [`schedule.ccl`](../../implementations/GRHayLib/schedule.ccl):
  `GRHayLib_initialize` scheduling at `CCTK_WRAGH` only when
  `ID_converter_ILGRMHD` is not active, and `GRHayLib_terminate` at
  `CCTK_TERMINATE`.
- [`src/GRHayLib.h`](../../implementations/GRHayLib/src/GRHayLib.h): public
  thorn header, aggregate includes, and external `ghl_params`/`ghl_eos`
  declarations.
- [`src/initialize_and_shutdown.c`](../../implementations/GRHayLib/src/initialize_and_shutdown.c):
  parameter validation, Con2Prim keyword parsing, EOS initialization,
  function-pointer selection, allocation, and cleanup.
- [`src/make.code.defn`](../../implementations/GRHayLib/src/make.code.defn):
  direct-compile source list and module subdirectory registry.

## Child Routes

- [Cactus build boundary](grhaylib/cactus-build-boundary.md): route for
  `configuration.ccl`, `interface.ccl`, `schedule.ccl`, `GRHayLib.h`,
  `src/make.code.defn`, HDF5, and direct-compile source layout.
- [Runtime parameter contract](grhaylib/runtime-parameter-contract.md):
  route for `GRHayLib_paramcheck`, `GRHayLib_initialize`,
  `GRHayLib_terminate`, `ghl_params`, `ghl_eos`, EOS setup, Con2Prim parser
  mapping, and Cactus parameter mapping.
- [Verification and drift](grhaylib/verification-and-drift.md): route
  for manual downstream verification, CI/test caveats, header aggregation
  checks, source-list checks, and parameter/parser drift.

## Gem Routes

Do not duplicate GRHayL physics or solver behavior here. Use GRHayLib only for
thorn ownership, Cactus parameter mapping, direct-compile build boundaries, and
global lifecycle.

- Core structs, shared parameters, EOS dispatch, errors, and utilities:
  [Core/Chalice](../core/index.md),
  [Shared parameters and enums](../core/shared-parameters-and-enums.md),
  [Core EOS dispatch](../core/eos-dispatch-contract.md).
- Atmosphere behavior: [Atmosphere hub](../gems/atmosphere.md).
- Con2Prim recovery and solver support:
  [Con2Prim recovery flow](../gems/con2prim/recovery-flow.md) and
  [Con2Prim solver matrix](../gems/con2prim/solver-matrix.md).
- EOS initialization, hybrid EOS, tabulated EOS, HDF5 table behavior:
  [EOS hub](../gems/eos.md),
  [EOS initialization and dispatch](../gems/eos/initialization-and-dispatch.md),
  [Tabulated EOS table contract](../gems/eos/tabulated-table-contract.md).
- Flux_Source behavior:
  [Flux Source hub](../gems/flux-source.md).
- Induction behavior:
  [Induction hub](../gems/induction.md).
- Neutrinos and NRPyLeakage behavior:
  [Neutrinos hub](../gems/neutrinos.md).
- Reconstruction behavior and PPM parameters:
  [Reconstruction hub](../gems/reconstruction.md) and
  [PPM flow](../gems/reconstruction/ppm-flow.md).

## Notes For Readers

- `GRHayLib.h` is the public thorn entry point, but upstream public API truth
  still lives in [`GRHayL/include/`](../../GRHayL/include/).
- `ghl_params` and `ghl_eos` are shared global thorn state. Treat lifecycle,
  pointer ownership, and tabulated EOS cleanup as GRHayLib runtime-contract
  concerns.
- Cactus parameter keywords and parser cases do not prove a Con2Prim method is
  fully supported. Compare against the Con2Prim solver matrix, selector
  dispatch, build lists, and tests.
- Current checkout proves static CCL/header/source-registry parity only. Copied
  module/include layout is absent; Cactus compile/runtime remain unverified.
