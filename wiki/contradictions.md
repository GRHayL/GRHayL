# Current Contradictions And Maintainer Decisions

Purpose: hold only unresolved current repository conflicts and bounded product
decisions. This is not history. Owner leaves document observed behavior and
coverage gaps; rows below name only disputed intent or unsafe seams needing a
product/docs choice. Delete a row after source-backed resolution.

Safe rule while any row remains: describe observable current behavior, never
upgrade declaration/source/workflow presence into built, initialized, executed,
or supported status.

| Area and bounded proposition | Competing evidence and observed behavior | Impact and safe pending wording | Smallest coherent decision |
| --- | --- | --- | --- |
| Configure build-type names | `configure` help advertises `nocflags`, which its parser rejects; parser accepts undocumented `plain`. | Document parser behavior as current and the no-flags help string as broken; users following it can fail configuration. | Choose the accepted no-flags name and align help and parser. |
| Legacy static Makefile route | [generate_makefile.sh](../generate_makefile.sh) scans all repo manifests, emits GRHayLib initialization as a target, and uses checkout-invalid `-I./include`. | Current output is broken; generation may be reproduced only in a disposable tree and `make` must not run. | Repair parser/scope/include root or retire the legacy route. |
| Fixture repository name | [README.md](../README.md) names `GRHayL/Test_Data`; [docs/raw/mainpage.md](../docs/raw/mainpage.md), workflows, and [.github/run_tests.sh](../.github/run_tests.sh) use operational `GRHayL/TestData`. | Use `TestData` for current scripted behavior; README link can misroute users. | Choose canonical repository name and align docs/scripts. |
| Perturbation formula | `README.md` states fixed `input * (1 + 1e-14)`; `docs/raw/mainpage.md` and `Unit_Tests/data_gen/*.c` use test-specific magnitudes and often random signed factors. | Say perturbation is generator-specific; fixed global formula is false. | Update README or standardize every generator. |
| CI path-filter policy | Every workflow ignores docs/wiki/Markdown for ignored-only push/PR changes. Ubuntu-GCC can select normal GRHayL jobs for implementation-only changes containing a non-ignored path, subject to branch/event filters, but does not inspect GRHayLib symlinks; the other workflows still ignore `implementations/**`. Schedules and mixed-path changes differ. `Makefile_old` is also ignored though absent. | Do not claim docs/KB CI proof or GRHayLib/Cactus validation; eligible implementation-only changes select ordinary GRHayL jobs; Markdown-only implementation changes remain ignored. | Decide whether broader implementation triggers or a real downstream build are wanted and whether the stale ignore entry remains intentional. |
| Coverage action meaning | Ubuntu clang jobs invoke coverage despite action note that expected files are absent; only some Intel jobs invoke it while Intel collection body is commented; macOS GCC invokes action with local body commented, macOS clang does not. Upload-action invocation alone does not prove that a usable coverage artifact was collected or uploaded. | Describe coverage as compiler/job-specific and unverified where collection is commented/missing. | Define intended artifact policy per compiler/OS, then align jobs and composite action. |
| macOS compiler setup | [Compile action](../.github/actions/compile_GRHayL/action.yml) has an unmatched quote inside the GCC `CC` command substitution; its clang branch supplies `--cflags` repeatedly, while [configure](../configure) keeps only the last value. | Workflow labels do not establish successful GCC selection; the later clang include flag replaces the requested coverage flags. | Correct GCC command quoting and pass the intended clang flags in one argument. |
| Codecov peer synchronization | [Ubuntu-GCC finalizer](../.github/workflows/github-actions-Ubuntu-gcc.yml) polls only returned, recent peer runs matching head/event and waits for `completed`, without requiring each named workflow to appear or checking conclusions. | Notification is not proof that every peer ran, succeeded, or uploaded coverage. | Decide whether completion of visible peers is sufficient or require explicit peer existence/success and upload evidence. |
| Radial atmosphere API | `ghl_set_prims_to_radial_falloff_atm` is installed-header declared and source-present, but absent from Atmosphere manifest/docs/tests; source leaves density/pressure/energy assignments commented. | No standard built runtime path; do not call it supported. | Complete/build/document/test routine or remove/internalize stale surface. |
| Palenzuela iteration diagnostics | [Public diagnostics contract](../GRHayL/include/ghl_con2prim.h) describes successful-solver counters, but [Brent](../GRHayL/Con2Prim/brent.c) returns early before initializing `n_iters`; shared Palenzuela callers copy it without prior initialization. | Successful recovery alone does not establish a valid count; a tabulated retry can retain a preceding attempt's count. See [internal numerics](gems/con2prim/internal-numerics-and-root-finding.md#brent-and-root-helper-route). | Initialize the counter for every root-finder call and align the documented diagnostic semantics. |
| Tabulated reinitialization in the error test | The public lifecycle requires cleanup before reinitializing a live table, while [Unit_Tests/unit_test_code_error.c](../Unit_Tests/unit_test_code_error.c) repeatedly initializes one live `tab_eos` object. | The cases can abandon table allocations; do not treat them as lifecycle examples or leak-check evidence. | Clean the live table before each reinitialization or give each case a fresh object. |
| Negative metric determinant | `ghl_initialize_metric` divides cofactors by `fabs(detgamma)`; determinant-enforced path preserves negative conformal determinant sign, warns, returns `void`, and test does not require positive determinant. | Validity requires positive-definite metric; negative inputs have no checked failure status and inverse semantics are unsafe. | Reject/report invalid metric or explicitly define recovery behavior and test invariant. |
| GRHayLib lifecycle handoff | Schedule skips initialize when `ID_converter_ILGRMHD` is active but always schedules terminate; globals use unchecked, uninitialized `malloc`; local repo contains no alternate owner or null/partial guards. | Normal success path is documented; skipped/failed initialization ownership is unverified and termination can dereference invalid state. | Define alternate owner/schedule dependency and add allocation/partial-state guards. |
| Reconstruction fixture provenance | ET_Legacy generator writes perturbed input that replay does not read; trusted and perturbed outputs remain externally supplied. | Generator existence does not prove use or provenance of every fixture path. | Consume or remove the unused perturbed-input path and document output provenance. |

## Owner Routes

- Build, runner, CI, installed surface: [Build And CI](build-and-ci.md),
  [Generated Boundaries](generated-boundaries.md), [Test Map](test-map.md).
- Core/Atmosphere: [Core](core/index.md),
  [Atmosphere prescription](gems/atmosphere/prescription-contract.md).
- Con2Prim/EOS: [solver matrix](gems/con2prim/solver-matrix.md),
  [EOS initialization](gems/eos/initialization-and-dispatch.md),
  [tabulated catalog](gems/eos/tabulated-interpolator-catalog.md).
- Flux/Induction/Reconstruction/Neutrinos: respective [gem router](gems/index.md).
- Downstream: [GRHayLib runtime contract](implementations/grhaylib/runtime-parameter-contract.md)
  and [verification boundary](implementations/grhaylib/verification-and-drift.md).
