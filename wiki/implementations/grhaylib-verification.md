# GRHayLib Verification

> Static parity checks, CI path-filter behavior, drift checklist, and external Cactus proof boundaries for `implementations/GRHayLib/`. · Status: confirmed
> Up: [Implementations](index.md)

## Summary

This repository supports static source, CCL, header-aggregation, symlink, and parameter/parser cross-checks for `implementations/GRHayLib/`, but contains no Cactus flesh or thorn build environment. Neither CI library jobs nor `Unit_Tests/` (`ET_Legacy` included) prove that GRHayLib compiles, schedules, or runs inside Cactus or the Einstein Toolkit.

## Detail

### CI Coverage And Proof Boundaries

Claim status: contested; contradiction: CONTR-0005. See [CONTR-0005](../contradictions.md#contr-0005).

For compiler-workflow `paths-ignore` selection on `implementations/**`, see [CI Workflows](../architecture/ci-workflows.md#compiler-and-os-workflow-matrix); those standalone library jobs and `.github/run_tests.sh` neither inspect GRHayLib symlinks nor compile the Cactus thorn.

| Claim | Verification boundary in this repository |
| --- | --- |
| CCL, `GRHayLib.h`, `SUBDIRS`, and parameter/parser parity | Checkable statically against `GRHayL/include/` and `GRHayL/**/make.code.defn`. (Ev: `implementations/GRHayLib/src/make.code.defn::SUBDIRS`) |
| Upstream source/header paths resolve under `src/` | Verified locally via tracked symlinks; does not prove Cactus build compatibility. (Ev: `implementations/GRHayLib/src/GRHayLib.h`) |
| Cactus compile, link, schedule, parameter steering, and runtime cleanup | Unverified here; requires an external Cactus/Einstein Toolkit checkout and build command. (Ev: `implementations/GRHayLib/schedule.ccl`) |
| Upstream `Unit_Tests/` or `ET_Legacy` prove GRHayLib | Inapplicable; they test standalone GRHayL C routines, not thorn integration. |

### Static Parity And Drift Checklist

Verify these items statically whenever `implementations/GRHayLib/` or shared `GRHayL/include/` declarations change:

1. **Header aggregation parity:** every header directly included by `src/GRHayLib.h` exists in `GRHayL/include/` and appears in `GRHayL/include/make.code.defn`, and the `src/include` symlink resolves.
2. **Direct-compile source-list parity:** every `SUBDIRS` entry in `src/make.code.defn` resolves via tracked symlinks to an upstream directory, and every upstream source-bearing manifest directory is listed.
3. **Parameter and parser parity:** compare `param.ccl` keywords and ranges against `parse_C2P_routine_keyword`, `parse_eos_table_type_keyword`, `GRHayLib_paramcheck`, and `GRHayLib_initialize` in `src/initialize_and_shutdown.c`.
4. **Upstream API drift points:** recheck `configuration.ccl` (`requires HDF5`), EOS initializer and `ghl_tabulated_free_memory` signatures, `ghl_parameters` and `ghl_eos_parameters` fields, `ghl_con2prim_id_t` support in [Solver Matrix](../gems/con2prim/solver-matrix.md), PPM overrides, and `CCTK_WRAGH`/`CCTK_TERMINATE` cleanup with `ID_converter_ILGRMHD`.
5. **Targeted upstream tests & docs guard:** run targeted `Unit_Tests/` only when upstream API behavior changed; compare `git status --short --untracked-files=all -- docs Doxyfile`, `git diff HEAD -- docs Doxyfile`, and `git diff --cached -- docs Doxyfile` with the task's starting state to confirm no unintended staged or unstaged changes and no new untracked doc files (pre-existing edits belong to the starting state). Never invent a local Cactus build command.

Evidence (descriptive behavior): `implementations/GRHayLib/src/initialize_and_shutdown.c::GRHayLib_paramcheck`; corroboration: `GRHayL/include/make.code.defn::install_headers`.

## Sources

- [implementations/GRHayLib/README](../../implementations/GRHayLib/README) - thorn scope
- [implementations/GRHayLib/configuration.ccl](../../implementations/GRHayLib/configuration.ccl) - `requires HDF5`
- [implementations/GRHayLib/interface.ccl](../../implementations/GRHayLib/interface.ccl) - thorn interface
- [implementations/GRHayLib/param.ccl](../../implementations/GRHayLib/param.ccl) - Cactus parameters
- [implementations/GRHayLib/schedule.ccl](../../implementations/GRHayLib/schedule.ccl) - `CCTK_WRAGH` / `CCTK_TERMINATE`
- [implementations/GRHayLib/src/GRHayLib.h](../../implementations/GRHayLib/src/GRHayLib.h) - aggregate header
- [implementations/GRHayLib/src/initialize_and_shutdown.c](../../implementations/GRHayLib/src/initialize_and_shutdown.c) - parameter validation, initialization, and cleanup code
- [implementations/GRHayLib/src/make.code.defn](../../implementations/GRHayLib/src/make.code.defn) - `SUBDIRS` manifest
- [GRHayL/include/make.code.defn](../../GRHayL/include/make.code.defn) - upstream `install_headers`

## See Also

- Depends on: [GRHayLib Cactus Build](grhaylib-cactus-build.md)
- Depends on: [GRHayLib Runtime Parameters](grhaylib-runtime-parameters.md)
- See also: [CI Workflows](../architecture/ci-workflows.md)
- See also: [Change Impact](../architecture/change-impact.md)
- See also: [Contradictions](../contradictions.md#contr-0005)
