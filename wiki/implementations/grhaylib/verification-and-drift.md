# GRHayLib Verification And Drift

This page routes verification for implementation-only changes under
`implementations/GRHayLib/` and upstream GRHayL changes that can drift the
GRHayLib Cactus thorn boundary.

Repo-local evidence is the authority here. Use this page to choose checks, not
to replace the implementation files or upstream tests.

## Coverage Boundary

Repo GitHub Actions ignore implementation-only pull requests and pushes:
workflow `paths-ignore` entries include `implementations/**`, along with
`docs/**` and `*.md`. The workflows run `ET-Legacy` jobs and other unit-test
jobs directly; `.github/run_tests.sh` is a separate scripted test driver. Those
jobs build and test the GRHayL library test programs, not a Cactus thorn.

`Unit_Tests/` does not directly cover GRHayLib/Cactus thorn behavior. Its tests
exercise upstream GRHayL source areas, fixture comparisons, generated data, and
ET_Legacy comparisons. ET_Legacy tests are related legacy comparisons for
upstream behavior; they are not proof that `implementations/GRHayLib/` builds,
schedules, initializes, parses Cactus parameters, or runs inside Cactus or the
Einstein Toolkit.

## Verification Hierarchy

1. Manual source cross-check. Compare claims against
   [README](../../../implementations/GRHayLib/README),
   [configuration.ccl](../../../implementations/GRHayLib/configuration.ccl),
   [interface.ccl](../../../implementations/GRHayLib/interface.ccl),
   [param.ccl](../../../implementations/GRHayLib/param.ccl),
   [schedule.ccl](../../../implementations/GRHayLib/schedule.ccl),
   [src/GRHayLib.h](../../../implementations/GRHayLib/src/GRHayLib.h),
   [src/initialize_and_shutdown.c](../../../implementations/GRHayLib/src/initialize_and_shutdown.c),
   and [src/make.code.defn](../../../implementations/GRHayLib/src/make.code.defn).
2. Manual repo-relative link check. Verify new wiki links resolve and avoid
   generated docs or external-only paths.
3. Header aggregation parity check. Current static comparison finds all ten
   headers directly included by `src/GRHayLib.h` in both `GRHayL/include/` and
   the upstream install manifest. Recheck after public-header changes. This
   does not create the absent `src/include/` copy layout.
4. Direct-compile source-list check. Current static comparison finds that all
   29 GRHayLib `SUBDIRS` name upstream directories and every one of the 28
   upstream source-bearing manifest directories is listed. The checkout still
   lacks those copied directories below the thorn, so this is registry parity,
   not a successful build or proof of an external copy process.
5. Parameter/parser parity check. Compare Cactus keywords in `param.ccl` with
   `parse_C2P_routine_keyword`, `parse_eos_table_type_keyword`,
   `GRHayLib_paramcheck`, and `GRHayLib_initialize` in
   `src/initialize_and_shutdown.c`.
6. Upstream GRHayL targeted tests. Run targeted `Unit_Tests/` only when an
   upstream public API or behavior change could affect GRHayLib. These tests
   can support upstream confidence; they do not verify Cactus thorn build or
   runtime behavior.
7. External/manual Cactus or Einstein Toolkit build. When a downstream
   Cactus/Einstein Toolkit checkout is available, perform its build and runtime
   verification there. Do not invent a local Cactus build command in this repo.
8. Generated-doc guard. `git diff -- docs` and `git diff -- Doxyfile` must be
   empty for implementation verification work unless a separate task explicitly
   asks for documentation-source changes.

## Drift Checklist

- HDF5: `configuration.ccl` has hard `requires HDF5`, while core GRHayL
  `configure` supports `--disable-hdf5`; do not assume no-HDF5 library checks
  prove GRHayLib thorn behavior.
- Headers: `src/GRHayLib.h` aggregates `ghl.h`, atmosphere, Con2Prim,
  induction, reconstruction, flux/source, radiation, hybrid/tabulated EOS, and
  NRPyLeakage headers. Header rename/split/addition can require parity review.
- Source list: `src/make.code.defn` must track upstream direct-compile source
  directories and new module subdirectories, with the absent-subdir caveat
  above.
- EOS signatures: check calls to `ghl_initialize_simple_eos_functions_and_params`,
  `ghl_initialize_hybrid_eos_functions_and_params`,
  `ghl_initialize_eos_functions`, `ghl_initialize_tabulated_eos`, and
  `ghl_tabulated_free_memory`.
- `ghl_parameters`: check `ghl_initialize_params`, Cactus parameter mapping,
  Con2Prim fields, entropy/temperature flags, Lorentz/gauge values, and PPM
  field copies.
- `ghl_eos_parameters`: check atmosphere, min/max, hybrid, tabulated,
  table-type, table-bounds, and cleanup expectations.
- `ghl_con2prim_id_t`: check enum values, parser keywords, backup routine
  slots, and invalid-key handling.
- Solver support: Cactus keywords and parser cases do not prove solver support;
  compare against Con2Prim dispatch, build-list evidence, tests, and
  `wiki/gems/con2prim/solver-matrix.md`.
- PPM fields: `param.ccl` exposes flattening and shock parameters that overwrite
  `ghl_initialize_params` defaults after initialization.
- Lifecycle: `GRHayLib_initialize` allocates global `ghl_params` and `ghl_eos`;
  `GRHayLib_terminate` frees tabulated memory first when needed, then frees
  both globals and clears them.
- Cactus schedule timing: initialization is scheduled at `CCTK_WRAGH` unless
  `ID_converter_ILGRMHD` is active; termination is scheduled at
  `CCTK_TERMINATE`. Global-pointer consumers depend on this timing.

## Current Proof State

| Claim | Current state |
| --- | --- |
| CCL/header/source-registry parity | Static source cross-check only; current names agree as bounded above. |
| Thorn files locally form a compilable copied layout | False for this checkout: module subdirectories and `src/include/` are absent. |
| Cactus compile/link | Unverified; no real Cactus build environment or owner command is supplied here. |
| Schedule/parameter/runtime behavior | Unverified; source describes intended Cactus behavior, but no Cactus execution evidence exists. |
| ET_Legacy or upstream Unit_Tests prove GRHayLib | False; they test upstream GRHayL behavior, not this thorn integration. |

## Local Ground Truth

- [Workflows](../../workflows.md)
- [Test Map](../../test-map.md)
- [Change Impact](../../change-impact.md)
- [Build And CI](../../build-and-ci.md)
- [.github/workflows](../../../.github/workflows/)
- [.github/run_tests.sh](../../../.github/run_tests.sh)
- [Unit_Tests](../../../Unit_Tests/)
- [implementations/GRHayLib/README](../../../implementations/GRHayLib/README)
- [implementations/GRHayLib/configuration.ccl](../../../implementations/GRHayLib/configuration.ccl)
- [implementations/GRHayLib/interface.ccl](../../../implementations/GRHayLib/interface.ccl)
- [implementations/GRHayLib/param.ccl](../../../implementations/GRHayLib/param.ccl)
- [implementations/GRHayLib/schedule.ccl](../../../implementations/GRHayLib/schedule.ccl)
- [implementations/GRHayLib/src/GRHayLib.h](../../../implementations/GRHayLib/src/GRHayLib.h)
- [implementations/GRHayLib/src/initialize_and_shutdown.c](../../../implementations/GRHayLib/src/initialize_and_shutdown.c)
- [implementations/GRHayLib/src/make.code.defn](../../../implementations/GRHayLib/src/make.code.defn)
