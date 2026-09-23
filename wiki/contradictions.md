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
| Legacy static Makefile route | `generate_makefile.sh` exits successfully but scans all repo manifests, emits GRHayLib initialization as a target, and uses checkout-invalid `-I./include`. | Current output is broken; generation may be reproduced only in a disposable tree and `make` must not run. | Repair parser/scope/include root or retire the legacy route. |
| Fixture repository name | `README.md` names `GRHayL/Test_Data`; `docs/raw/mainpage.md`, workflows, and `.github/run_tests.sh` use operational `GRHayL/TestData`. | Use `TestData` for current scripted behavior; README link can misroute users. | Choose canonical repository name and align docs/scripts. |
| Perturbation formula | `README.md` states fixed `input * (1 + 1e-14)`; `docs/raw/mainpage.md` and `Unit_Tests/data_gen/*.c` use test-specific magnitudes and often random signed factors. | Say perturbation is generator-specific; fixed global formula is false. | Update README or standardize every generator. |
| CI path-filter policy | Every workflow ignores docs/wiki/Markdown for ignored-only push/PR changes. Ubuntu-GCC runs normal GRHayL jobs for implementation-only changes but does not inspect GRHayLib symlinks; the other four workflows still ignore `implementations/**`. Schedules and mixed-path changes differ. `Makefile_old` is also ignored though absent. | Do not claim docs/KB CI proof or GRHayLib/Cactus validation; implementation-only changes receive only ordinary GRHayL job evidence. | Decide whether broader implementation triggers or a real downstream build are wanted and whether the stale ignore entry remains intentional. |
| Coverage action meaning | Ubuntu clang jobs invoke coverage despite action note that expected files are absent; only 2/13 Intel jobs invoke it while Intel collection body is commented; macOS GCC invokes action with local body commented, macOS clang does not. Upload-action invocation alone does not prove that a usable coverage artifact was collected or uploaded. | Describe coverage as compiler/job-specific and unverified where collection is commented/missing. | Define intended artifact policy per compiler/OS, then align jobs and composite action. |
| Radial atmosphere API | `ghl_set_prims_to_radial_falloff_atm` is installed-header declared and source-present, but absent from Atmosphere manifest/docs/tests; source leaves density/pressure/energy assignments commented. | No standard built runtime path; do not call it supported. | Complete/build/document/test routine or remove/internalize stale surface. |
| `tau_atm` family semantics | Simple/hybrid initialize `rho_atm * eps_atm`; tabulated initializes `rho_min * eps_min`. Public struct calls it an atmosphere value; Con2Prim uses it as floor. | Formulas are not assumed equivalent. | Decide whether field is atmosphere energy or minimum conservative floor and align names/formulas/docs. |
| Tabulated reinitialization in the error test | The public lifecycle requires cleanup before reinitializing a live table, while `Unit_Tests/unit_test_code_error.c` repeatedly initializes one live `tab_eos` object. | The cases can abandon table allocations; do not treat them as lifecycle examples or leak-check evidence. | Clean the live table before each reinitialization or give each case a fresh object. |
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
