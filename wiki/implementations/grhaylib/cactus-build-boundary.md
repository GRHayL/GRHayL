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

All statements above are static CCL observations. This checkout supplies no
Cactus flesh/configuration, thorn-list build, parameter file, or runtime result.
ET_Legacy unit tests are upstream regression comparisons and cannot substitute
for a real Cactus build or scheduled run.

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

Static parity in this checkout: all ten headers named directly by `GRHayLib.h`
exist under `GRHayL/include/` and all ten appear in the upstream install-header
manifest. `ghl.h` also supplies several transitive Core includes. This confirms
names only; the thorn header spells them as `./include/...`, and no such copied
`implementations/GRHayLib/src/include/` tree exists locally.

## Thorn Source Registry

`implementations/GRHayLib/src/make.code.defn` is the GRHayLib thorn source
registry/direct-compile boundary. It lists local `SRCS` as:

- `initialize_and_shutdown.c`

It also lists `SUBDIRS` for module paths including `Atmosphere`, `Con2Prim`,
`Con2Prim/Hybrid`, `Con2Prim/Tabulated`, `EOS/Hybrid`, `EOS/Tabulated`,
`Flux_Source`, `GRHayL_Core`, `Induction`, `Neutrinos/NRPyLeakage`, and
`Reconstruction` variants.

Static registry parity is complete at directory level: each of the 29 listed
`SUBDIRS` exists under upstream `GRHayL/`, and every one of the 28 upstream
directories whose reachable manifest has a nonempty `SRCS` block is listed.
`Con2Prim/Hybrid` is the one extra intermediate directory. This proves current
registry-name agreement, not that Cactus received copied source.

Absent-copy-layout caveat: this checkout has only `doc/` and `src/` below
`implementations/GRHayLib/`. None of the 29 listed module subdirectories and no
`src/include/` aggregate-header dependency are present there. Therefore the
registry is a direct-compile/copy-layout contract, not a locally complete thorn
build tree. Do not call the absent layout a successful build, describe it as
generated output, or infer a copy step not present in repository source.

## Verification Status

- Static CCL/header/source-registry parity: checked against this checkout.
- Local thorn layout: incomplete, so no local Cactus compile claim.
- Cactus build, schedule execution, parameter parsing, and runtime cleanup:
  unverified; require a real Cactus/Einstein Toolkit environment and its
  owner-provided command.
- Core GRHayL and ET_Legacy test results: inapplicable as proof of those thorn
  states.

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
