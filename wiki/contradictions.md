# Contradictions And Drift

Repo ground truth: `README.md`, `docs/raw/mainpage.md`,
`.github/actions/code-coverage/action.yml`, `.github/workflows/`,
`.github/run_tests.sh`, and `Unit_Tests/`.

This page records open drift points and questions. It is not a running history.

| Status | Area | Drift point | Evidence | Next action |
| --- | --- | --- | --- | --- |
| `open` | Test data repository name | `README.md` links `GRHayL/Test_Data`; `docs/raw/mainpage.md`, `.github/run_tests.sh`, and workflows use `GRHayL/TestData`. | `README.md` CI Testing section; `docs/raw/mainpage.md` CI Testing section; `.github/run_tests.sh` `test_data_base_url`. | Choose one canonical fixture repository name and update docs/scripts if needed. |
| `open` | Perturbation wording | `README.md` says perturbed input is `input X (1 + 1e-14)`. `docs/raw/mainpage.md` says a small test-specific relative factor, commonly with `rand(-1,1) * epsilon`. Data generators show multiple perturbation magnitudes and random factors. | `README.md` CI Testing section; `docs/raw/mainpage.md` CI Testing section; `Unit_Tests/data_gen/*.c`. | Prefer docs wording that reflects generator behavior; avoid a single fixed perturbation rule unless all generators match it. |
| `open` | Configure build types | `configure -h` lists `nocflags`, but the implemented build-type case accepts `plain` and does not accept `nocflags`. | `configure` help text and build-type `case` block. | Confirm intended no-flags build type and align help or parser. |
| `open` | Coverage expectations for Linux clang | Coverage action says Linux clang still does not generate expected coverage files, but Ubuntu clang workflows call the coverage action. | `.github/actions/code-coverage/action.yml`; `.github/workflows/github-actions-Ubuntu-clang.yml`. | Confirm whether clang coverage upload is expected to be useful or only best-effort. |
| `needs maintainer` | Coverage expectations for Linux Intel | Intel coverage body is commented in the action, but some Ubuntu Intel jobs still invoke the coverage action. | `.github/actions/code-coverage/action.yml`; `.github/workflows/github-actions-Ubuntu-intel.yml`. | Decide whether Intel coverage jobs should be disabled, implemented, or documented as upload-only/no local report. |
| `needs maintainer` | Coverage expectations for macOS | macOS coverage command bodies are commented with notes about tool selection, while macOS GCC workflow invokes coverage and macOS clang has coverage steps commented out. | `.github/actions/code-coverage/action.yml`; `.github/workflows/github-actions-MacOS-gcc.yml`; `.github/workflows/github-actions-MacOS-clang.yml`. | Confirm intended macOS coverage policy. |
| `open` | Docs-only CI behavior | Workflows ignore `docs/**` and `*.md` on push/pull_request, so docs-only changes may not receive CI validation. | `paths-ignore` in every workflow. | Decide whether docs checks should run elsewhere or remain manual. |
| `needs maintainer` | Core test fixture spelling | `Unit_Tests/data_gen/unit_test_data_grhayl_core_test_suite.c` writes `grhayL_core_test_suite_input.bin`, while scripts/tests use `grhayl_core_test_suite_input.bin`. | Data generator string; `.github/run_tests.sh`; `Unit_Tests/unit_test_grhayl_core_test_suite.c`. | Confirm whether generator output name is stale or intentionally unused. |

## Resolved

No resolved drift points recorded in this page yet.
