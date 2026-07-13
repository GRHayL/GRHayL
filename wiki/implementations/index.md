# Implementation KB Router

This page routes implementation KB work for `implementations/*`. It covers
downstream trees that interface with GRHayL by direct-compile integration, not
the normal installed-library path described by the root README.

Repo-local implementation files remain the ground truth. KB pages here point
readers to implementation sources and Cactus metadata; they do not replace
`implementations/*` files or upstream GRHayL source.

## Scope

- `implementations/*` is the implementation scope for this router.
- [GRHayLib](grhaylib.md) is currently the only implementation tree present:
  [`implementations/GRHayLib/`](../../implementations/GRHayLib/).
- Direct-compile implementation paths are for infrastructures that compile
  GRHayL code inside their own build system. Normal users usually link against
  a configured and installed GRHayL library instead.
- KB routing does not imply changes to `docs/**`, generated docs, `Doxyfile`,
  or `implementations/**`.

## Routes

| Page | Use it for | Ground truth |
| --- | --- | --- |
| [GRHayLib hub](grhaylib.md) | Overview for the Cactus/Einstein Toolkit thorn boundary. | [`implementations/GRHayLib/`](../../implementations/GRHayLib/) |
| [Cactus build boundary](grhaylib/cactus-build-boundary.md) | CCL metadata, exported header, direct-compile source list, HDF5 requirement, and source-layout caveats. | [`configuration.ccl`](../../implementations/GRHayLib/configuration.ccl), [`interface.ccl`](../../implementations/GRHayLib/interface.ccl), [`schedule.ccl`](../../implementations/GRHayLib/schedule.ccl), [`src/make.code.defn`](../../implementations/GRHayLib/src/make.code.defn) |
| [Runtime parameter contract](grhaylib/runtime-parameter-contract.md) | Global state, lifecycle, EOS setup, Con2Prim routing, and Cactus parameter mapping. | [`param.ccl`](../../implementations/GRHayLib/param.ccl), [`src/initialize_and_shutdown.c`](../../implementations/GRHayLib/src/initialize_and_shutdown.c), [`src/GRHayLib.h`](../../implementations/GRHayLib/src/GRHayLib.h) |
| [Verification and drift](grhaylib/verification-and-drift.md) | Manual downstream verification and implementation drift checks. | [`README`](../../implementations/GRHayLib/README), [`doc/documentation.tex`](../../implementations/GRHayLib/doc/documentation.tex), repo CI/test maps |

Current boundary: source registry and aggregate-header names have static parity
with upstream files, but copied thorn module/include layout is absent locally.
No Cactus build or runtime result is established by this repository's
ET_Legacy/upstream tests.

## Source Notes

- [`README.md`](../../README.md) distinguishes installed-library use from
  implementation trees that directly compile GRHayL for infrastructures such as
  the Einstein Toolkit.
