# Contradictions

> Register of active contested and stale GRHayL KB claims. · Status: confirmed

Active source discrepancies and unreconciled drift are tracked below. Recording an entry documents current code versus documentation or workflow behavior without implying a source patch.

## Rules

- Follow [KB Rules](SCHEMA.md#contradiction-rules): keep only `contested` and `stale` claims in `## Register`, paired with `### CONTR-XXXX` subsections and literal inline markers (`Claim status: contested; contradiction: CONTR-0001.` or `stale`) on every affected leaf.
- Look up one entry with `rg -n -F '| CONTR-0001 |' wiki/contradictions.md` and read its `### CONTR-XXXX` subsection.
- Retire resolved entries in the resolving commit once affected leaves are reconciled and all `#contr-*` links are removed. Retired IDs: none.

## Register

| ID | Claim status | Affected pages | Summary |
| --- | --- | --- | --- |
| CONTR-0001 | contested | [Build And Configure](architecture/build-and-configure.md) | `configure` help advertises `nocflags` (rejected); parser accepts undocumented `plain`. |
| CONTR-0002 | stale | [Build](architecture/build-and-configure.md); [Generated](architecture/generated-boundaries.md) | `generate_makefile.sh` scans outside `GRHayL/`, includes GRHayLib, emits `-I./include`. |
| CONTR-0003 | stale | [Ref](tests/reference-data-and-test-drivers.md), [Run](tests/runner-and-build-outputs.md), [ET](tests/et-legacy-comparison.md), [Cat](tests/unit-test-catalog.md) | `Test_Data` vs `TestData`. |
| CONTR-0004 | contested | [RefData](tests/reference-data-and-test-drivers.md); [ET](tests/et-legacy-comparison.md) | `README.md` states fixed `input * (1 + 1e-14)`; `data_gen/*.c` use test-specific random scales. |
| CONTR-0005 | contested | [CI](architecture/ci-workflows.md); [GRHayLib](implementations/grhaylib-verification.md) | Ubuntu-GCC `paths-ignore` omits `implementations/**` unlike peers; ignores absent `Makefile_old`. |
| CONTR-0006 | stale | [CI Workflows](architecture/ci-workflows.md) | Coverage action runs on clang, Intel, and macOS-GCC jobs where collection is commented or files are absent. |
| CONTR-0007 | contested | [CI Workflows](architecture/ci-workflows.md) | `compile_GRHayL/action.yml` has unmatched macOS GCC quote and repeated clang `--cflags` overwritten by `configure`. |
| CONTR-0008 | contested | [CI Workflows](architecture/ci-workflows.md) | Ubuntu-GCC `codecov-finalize` polls recent peer runs without checking peer presence, conclusions, or uploads. |
| CONTR-0009 | stale | [Atmosphere](gems/atmosphere.md); [API](architecture/public-api-map.md) | `ghl_set_prims_to_radial_falloff_atm` is declared with incomplete source, omitted from manifest/tests. |
| CONTR-0010 | contested | [TabLayout](gems/eos/tabulated-table-layout.md); [Errors](tests/expected-failure-and-error-keys.md) | `unit_test_code_error.c` reinitializes one live `tab_eos` without freeing table memory. |
| CONTR-0011 | contested | [Metric ADM](core/metric-adm.md) | `ghl_initialize_metric` divides by `fabs(detgamma)`; determinant-enforced init warns (`void`) on negative det. |
| CONTR-0012 | stale | [Recon](gems/reconstruction/tests-and-reference-data.md); [ET](tests/et-legacy-comparison.md) | ET_Legacy reconstruction generator writes perturbed input unused by replay; outputs are external. |

### CONTR-0001

- Competing sources & authority: `configure` help text advertising `--buildtype=nocflags` versus the `configure` argument parser rejecting `nocflags` and accepting undocumented `plain`. The executed parser decides current behavior; [Build And Configure](architecture/build-and-configure.md) is `contested`.
- Recheck trigger: edits to `configure` `--help` output or `--buildtype` `case` handling.
- Resolution check: run `./configure --help` and inspect the `--buildtype` `case` statement, verifying that the documented no-flags option matches the accepted parser token.

Evidence (descriptive behavior): `configure::nocflags`; corroboration: `configure::plain`.

### CONTR-0002

- Competing sources & authority: `generate_makefile.sh` scanning every repository `make.code.defn` (including `implementations/GRHayLib/src/initialize_and_shutdown.c`) and emitting `-I./include` versus the supported `configure` build path under `GRHayL/include/`. `configure` decides supported builds; [Build And Configure](architecture/build-and-configure.md) is `contested` while [Generated Boundaries](architecture/generated-boundaries.md) stays `confirmed`.
- Recheck trigger: edits to `generate_makefile.sh` or its removal from the repository.
- Resolution check: inspect `generate_makefile.sh` (or confirm removal) in a disposable checkout, verifying `GRHayL/`-scoped manifest discovery, `-IGRHayL/include` include root, shell shebang, and HDF5 filtering.

Evidence (descriptive behavior): `generate_makefile.sh::libgrhayl.a`; corroboration: `configure::GRHayL/make.code.defn`.

### CONTR-0003

- Competing sources & authority: `README.md` naming `GRHayL/Test_Data` versus `docs/raw/mainpage.md`, `.github/run_tests.sh`, and `.github/workflows/github-actions-Ubuntu-gcc.yml` using `GRHayL/TestData`. Operational scripts decide current download paths; [Reference Data And Test Drivers](tests/reference-data-and-test-drivers.md) and [ET Legacy Comparison](tests/et-legacy-comparison.md) are `contested` while [Runner And Build Outputs](tests/runner-and-build-outputs.md) and [Unit Test Catalog](tests/unit-test-catalog.md) stay `confirmed`.
- Recheck trigger: edits to `README.md`, `docs/raw/mainpage.md`, `.github/run_tests.sh`, or workflow TestData URLs.
- Resolution check: run `rg -n 'Test_Data|TestData' README.md docs/raw/mainpage.md .github/` and confirm one canonical repository name across docs and scripts.

Evidence (descriptive behavior): `README.md::Test_Data`; corroboration: `.github/run_tests.sh::TestData` and `docs/raw/mainpage.md::TestData`.

### CONTR-0004

- Competing sources & authority: `README.md` stating a uniform `input * (1 + 1e-14)` perturbation versus `docs/raw/mainpage.md` and `Unit_Tests/data_gen/*.c` using generator-specific magnitudes and random signed factors. Generator C sources decide actual perturbation math; affected leaves are `contested`.
- Recheck trigger: edits to `README.md` perturbation prose or `Unit_Tests/data_gen/unit_test_data_*.c`.
- Resolution check: compare `README.md` against `docs/raw/mainpage.md` and `Unit_Tests/data_gen/*.c`, confirming accurate description of generator-specific perturbations.

Evidence (descriptive behavior): `README.md::1e-14`; corroboration: `Unit_Tests/data_gen/unit_test_data_hybrid_flux.c`.

### CONTR-0005

- Competing sources & authority: `.github/workflows/github-actions-Ubuntu-gcc.yml` omitting `implementations/**` from `paths-ignore` (and ignoring absent `Makefile_old`) versus the peer compiler workflows ignoring `implementations/**`, while no workflow builds or tests the Cactus thorn. Workflow YAML decides selection; [CI Workflows](architecture/ci-workflows.md) is `contested` while [GRHayLib Verification](implementations/grhaylib-verification.md) stays `confirmed`.
- Recheck trigger: edits to `paths-ignore` in `.github/workflows/github-actions-*.yml`.
- Resolution check: inspect `paths-ignore` across every `.github/workflows/github-actions-*.yml` file and verify consistent `implementations/**` and `Makefile_old` policy.

Evidence (CI behavior): `.github/workflows/github-actions-Ubuntu-gcc.yml::Makefile_old`; corroboration: `.github/workflows/github-actions-Ubuntu-clang.yml::implementations/**`.

### CONTR-0006

- Competing sources & authority: `.github/workflows/github-actions-*.yml` invoking `.github/actions/code-coverage/action.yml` on Ubuntu clang, selected Ubuntu Intel, and macOS GCC jobs versus commented-out collection bodies and missing-file notes inside `code-coverage/action.yml`. Action shell steps decide actual collection; [CI Workflows](architecture/ci-workflows.md) is `contested`.
- Recheck trigger: edits to `.github/actions/code-coverage/action.yml` or coverage job steps in `.github/workflows/github-actions-*.yml`.
- Resolution check: inspect `code-coverage/action.yml` and every workflow job invoking it, verifying active collection commands or removal of no-op invocations.

Evidence (CI behavior): `.github/actions/code-coverage/action.yml::codecov/codecov-action@v5`; corroboration: `.github/workflows/github-actions-MacOS-gcc.yml`.

### CONTR-0007

- Competing sources & authority: `.github/actions/compile_GRHayL/action.yml` macOS GCC `CC` command substitution (unmatched quote in `egrep`) and repeated macOS clang `--cflags` arguments versus `configure` retaining only the final `--cflags` value. `configure` and shell syntax decide actual execution; [CI Workflows](architecture/ci-workflows.md) is `contested`.
- Recheck trigger: edits to `.github/actions/compile_GRHayL/action.yml` or `configure` `--cflags` parsing.
- Resolution check: inspect macOS GCC quoting and clang `--cflags` consolidation in `.github/actions/compile_GRHayL/action.yml` and confirm `configure` preserves all intended flags.

Evidence (CI behavior): `.github/actions/compile_GRHayL/action.yml::--cflags`; corroboration: `configure::--cflags`.

### CONTR-0008

- Competing sources & authority: `codecov.yml` deferred notifications and Ubuntu-GCC `codecov-finalize` polling returned recent peer runs for `completed` without requiring every named peer workflow to appear, succeed, or upload coverage. Workflow script decides actual synchronization; [CI Workflows](architecture/ci-workflows.md) is `contested`.
- Recheck trigger: edits to `codecov-finalize` in `.github/workflows/github-actions-Ubuntu-gcc.yml` or `codecov.yml`.
- Resolution check: inspect `codecov-finalize` polling logic and verify whether missing peers, failed conclusions, and upload confirmation are explicitly handled.

Evidence (CI behavior): `.github/workflows/github-actions-Ubuntu-gcc.yml::codecov-finalize`; corroboration: `codecov.yml`.

### CONTR-0009

- Competing sources & authority: `GRHayL/include/ghl_atmosphere.h` declaring `ghl_set_prims_to_radial_falloff_atm` and `GRHayL/Atmosphere/set_prims_to_radial_falloff_atm.c` leaving density/pressure/energy assignments commented versus `GRHayL/Atmosphere/make.code.defn` omitting that C file from built sources. Build manifest and C body decide current runtime availability; [Atmosphere](gems/atmosphere.md) is `stale` while [Public API Map](architecture/public-api-map.md) stays `confirmed`.
- Recheck trigger: edits to `ghl_atmosphere.h`, `set_prims_to_radial_falloff_atm.c`, or `GRHayL/Atmosphere/make.code.defn`.
- Resolution check: confirm either that `set_prims_to_radial_falloff_atm.c` is completed, listed in `make.code.defn`, documented, and tested, or that the unfinished declaration and file are removed.

Evidence (descriptive behavior): `GRHayL/include/ghl_atmosphere.h::ghl_set_prims_to_radial_falloff_atm`; corroboration: `GRHayL/Atmosphere/make.code.defn` and `GRHayL/Atmosphere/set_prims_to_radial_falloff_atm.c`.

### CONTR-0010

- Competing sources & authority: `ghl_tabulated_free_memory` cleanup requirements before reinitializing a live `ghl_eos_parameters` tabulated object versus `Unit_Tests/unit_test_code_error.c` repeatedly calling tabulated initialization on one live `tab_eos` object without freeing memory between error cases. Public EOS memory rules decide safe caller usage; affected leaves stay `confirmed`.
- Recheck trigger: edits to `Unit_Tests/unit_test_code_error.c` or tabulated EOS initialization cleanup.
- Resolution check: inspect `Unit_Tests/unit_test_code_error.c` and verify every live `tab_eos` object is freed before reinitialization or isolated per error case.

Evidence (descriptive behavior): `Unit_Tests/unit_test_code_error.c::tab_eos`; corroboration: `GRHayL/include/ghl_nrpyeos_tabulated.h::NRPyEOS_free_memory`.

### CONTR-0011

- Competing sources & authority: physical requirement for a positive-definite 3-metric versus `ghl_initialize_metric` dividing cofactors by `fabs(detgamma)` and `ghl_enforce_detgtij_and_initialize_ADM_metric` issuing a `ghl_warn` while returning `void` and preserving negative determinant sign. C implementations decide current behavior; [Metric And ADM Helpers](core/metric-adm.md) is `contested`.
- Recheck trigger: edits to `GRHayL/GRHayL_Core/initialize_metric.c`, `enforce_detgij_and_initialize_ADM_metric.c`, or `Unit_Tests/unit_test_grhayl_core_test_suite.c`.
- Resolution check: inspect both metric initializers and verify either checked error returns on nonpositive determinants or an explicit documented and tested recovery rule.

Evidence (descriptive behavior): `GRHayL/GRHayL_Core/initialize_metric.c::ghl_initialize_metric`; corroboration: `GRHayL/GRHayL_Core/enforce_detgij_and_initialize_ADM_metric.c::ghl_enforce_detgtij_and_initialize_ADM_metric`.

### CONTR-0012

- Competing sources & authority: `Unit_Tests/data_gen/unit_test_data_ET_Legacy_reconstruction.c` writing perturbed input data versus `Unit_Tests/unit_test_ET_Legacy_reconstruction.c` never reading that perturbed input while relying on externally supplied trusted and perturbed outputs. Test C files decide actual replay behavior; [ET Legacy Comparison](tests/et-legacy-comparison.md) is `contested` while [Reconstruction Tests And Reference Data](gems/reconstruction/tests-and-reference-data.md) stays `confirmed`.
- Recheck trigger: edits to `Unit_Tests/data_gen/unit_test_data_ET_Legacy_reconstruction.c` or `Unit_Tests/unit_test_ET_Legacy_reconstruction.c`.
- Resolution check: inspect both files and verify either that replay consumes the generated perturbed input or that unused generation code is removed and output provenance is documented.

Evidence (descriptive behavior): `Unit_Tests/unit_test_ET_Legacy_reconstruction.c`; corroboration: `Unit_Tests/data_gen/unit_test_data_ET_Legacy_reconstruction.c`.
