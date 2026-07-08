# GRHayLib Cactus Build Boundary

Purpose: route Cactus CCL, thorn header, and direct-compile source-list
questions for `implementations/GRHayLib/`. Use this page for GRHayLib thorn
build/interface boundaries; route upstream library build behavior through
[Build And CI](../../build-and-ci.md), generated-output questions through
[Generated Boundaries](../../generated-boundaries.md), and public header impact
through [Public API Map](../../public-api-map.md).

Sibling routes: [GRHayLib overview](../grhaylib.md),
[runtime parameter contract](runtime-parameter-contract.md), and
[verification and drift](verification-and-drift.md).

## CCL Boundary

`implementations/GRHayLib/configuration.ccl` is the Cactus configuration
boundary for this thorn. Local source states:

- `PROVIDES GRHayL`
- `requires HDF5`

That `requires HDF5` line is a hard GRHayLib thorn requirement in this local
implementation. Keep it separate from core GRHayL configured builds: root
`README.md`, `docs/raw/mainpage.md`, and `wiki/build-and-ci.md` describe core
GRHayL `--disable-hdf5` support through `configure`, with source/test filtering
and `GHL_DISABLE_HDF5`. The Cactus thorn file does not expose that no-HDF5
configuration path.

`implementations/GRHayLib/interface.ccl` is the Cactus interface boundary.
Local source states:

- `implements: GRHayLib`
- `INCLUDE HEADER: GRHayLib.h in GRHayLib.h`

Treat the included `GRHayLib.h` as the public Cactus-facing aggregate header,
not as a replacement for upstream module ownership in `GRHayL/include/`.

`implementations/GRHayLib/schedule.ccl` is the Cactus lifecycle boundary.
Local source schedules:

- `GRHayLib_initialize` at `CCTK_WRAGH`
- only when `ID_converter_ILGRMHD` is not active
- `GRHayLib_terminate` at `CCTK_TERMINATE`

These schedule names define thorn startup/shutdown placement for downstream
Cactus users. They do not prove standalone repo CI exercises GRHayLib; use
[verification and drift](verification-and-drift.md) for that boundary.

## Aggregate Header Boundary

`implementations/GRHayLib/src/GRHayLib.h` aggregates public GRHayL headers for
the thorn and declares the exported global pointers:

- `ghl_eos_parameters *ghl_eos`
- `ghl_parameters *ghl_params`

Header rename, split, include-order, struct, enum, or global-lifecycle changes
can become downstream API risk because Cactus thorns can include this aggregate
header rather than individual upstream headers. Route ownership of the upstream
headers through [Public API Map](../../public-api-map.md); route lifecycle and
parameter details through [runtime parameter contract](runtime-parameter-contract.md).

## Thorn Source Registry

`implementations/GRHayLib/src/make.code.defn` is the GRHayLib thorn source
registry/direct-compile boundary. It lists local `SRCS` as:

- `initialize_and_shutdown.c`

It also lists `SUBDIRS` for module paths including `Atmosphere`, `Con2Prim`,
`Con2Prim/Hybrid`, `Con2Prim/Tabulated`, `EOS/Hybrid`, `EOS/Tabulated`,
`Flux_Source`, `GRHayL_Core`, `Induction`, `Neutrinos/NRPyLeakage`, and
`Reconstruction` variants.

Absent-subdir caveat: this checkout has only `doc/` and `src/` below
`implementations/GRHayLib/` when checked with
`find implementations/GRHayLib -maxdepth 2 -type d | sort`. Therefore the
`SUBDIRS` list must be treated as a direct-compile/copy-layout boundary needing
maintainer confirmation, not proof those subdirectories exist locally. Do not
describe those listed paths as generated output, and do not infer a local
Cactus copy step unless source evidence is added.

## Local Ground Truth

- [`wiki/build-and-ci.md`](../../build-and-ci.md)
- [`wiki/generated-boundaries.md`](../../generated-boundaries.md)
- [`wiki/public-api-map.md`](../../public-api-map.md)
- [`implementations/GRHayLib/configuration.ccl`](../../../implementations/GRHayLib/configuration.ccl)
- [`implementations/GRHayLib/interface.ccl`](../../../implementations/GRHayLib/interface.ccl)
- [`implementations/GRHayLib/schedule.ccl`](../../../implementations/GRHayLib/schedule.ccl)
- [`implementations/GRHayLib/src/GRHayLib.h`](../../../implementations/GRHayLib/src/GRHayLib.h)
- [`implementations/GRHayLib/src/make.code.defn`](../../../implementations/GRHayLib/src/make.code.defn)
- [`README.md`](../../../README.md)
- [`docs/raw/mainpage.md`](../../../docs/raw/mainpage.md)
