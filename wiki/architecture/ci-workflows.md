# CI Workflows

> GitHub Actions compiler/OS matrices, job routing, composite actions, Codecov configuration, and `.github/run_tests.sh` overview. · Status: contested
> Up: [Architecture](index.md)

## Summary

GRHayL defines compiler/OS GitHub Actions workflows under `.github/workflows/`, composite actions under `.github/actions/`, root `codecov.yml`, and the scripted local replay driver `.github/run_tests.sh`. Workflow YAML establishes job selection only; several coverage bodies are commented out, macOS compiler setup has known quoting/flag defects, and Ubuntu-GCC `codecov-finalize` polls without checking peer presence or success.

## Detail

### Compiler And OS Workflow Matrix

| Workflow | Compiler | OS matrix | Coverage invocation |
| --- | --- | --- | --- |
| `github-actions-Ubuntu-gcc.yml` | `gcc` | `ubuntu-22.04`, `ubuntu-24.04` | Existing job groups invoke `code-coverage`; `compose-regularized-eos` uploads `compose-coverage.xml` directly. (Ev: `github-actions-Ubuntu-gcc.yml::compose-regularized-eos`) |
| `github-actions-Ubuntu-clang.yml` | `clang` | `ubuntu-22.04`, `ubuntu-24.04` | Every job invokes `code-coverage`; `c2p-failure` Ubuntu 24.04 `c2p_nn_guess` runs `--disable-hdf5` and inspects tabulated HLLE symbols without running a tabulated flux test. (Ev: `github-actions-Ubuntu-clang.yml::c2p-failure`) |
| `github-actions-Ubuntu-intel.yml` | `intel` / `icx` | `ubuntu-22.04`, `ubuntu-24.04` | Selected jobs invoke `code-coverage`; Intel collection body is commented out. (Ev: `github-actions-Ubuntu-intel.yml`) |
| `github-actions-MacOS-gcc.yml` | Homebrew GCC | `macos-15`, `macos-26` | Jobs invoke `code-coverage`; local collection body is commented out. (Ev: `github-actions-MacOS-gcc.yml`) |
| `github-actions-MacOS-clang.yml` | Homebrew LLVM clang | `macos-15`, `macos-26` | No job invokes `code-coverage`. (Ev: `github-actions-MacOS-clang.yml`) |

Claim status: contested; contradiction: CONTR-0005. See [CONTR-0005](../contradictions.md#contr-0005).

Every compiler workflow restricts `push` to `main`, has no `pull_request` branch filter, runs on cron `33 15 1,15 * *`, and skips push/PR runs only when **all** changed paths match `paths-ignore` (`docs/**`, `wiki/**`, Markdown/reStructuredText patterns, and absent `Makefile_old`). Ubuntu-GCC no longer ignores `implementations/**`, so eligible non-Markdown changes there (such as `implementations/GRHayLib/README` or `doc/documentation.tex`) select normal core library jobs without inspecting GRHayLib symlinks or building Cactus; peer workflows still ignore `implementations/**`. Mixed changes containing any non-ignored path and scheduled runs bypass `paths-ignore`.

Evidence (CI behavior): `github-actions-Ubuntu-gcc.yml::paths-ignore`; corroboration: `github-actions-Ubuntu-clang.yml::implementations/**`.

### Job Routing Across Test Groups

| Job | Scoped tests |
| --- | --- |
| `ET-Legacy` | `conservs`, `primitives`, `induction_gauge_rhs`, `HLL_flux`, `reconstruction`, `flux_source` (Ev: `github-actions-Ubuntu-gcc.yml::ET-Legacy`) |
| `c2p-routines` / `c2p-failure` | Conservative limits, hybrid multi-method, primitive limits/`u0`, conservs/`Tmunu`; `hybrid_failure`, `c2p_nn_guess` |
| `tabulated-eos` / `piecewise-polytrope-eos` / `con-to-prim-tabulated` | Tabulated EOS read/interpolation, piecewise-polytrope EOS, tabulated Con2Prim |
| `grhayl-core` / `code-failure` | Core EOS/default, atmosphere, metric, routine-name suite; expected error-code keys |
| `flux` / `reconstruction` | `hybrid_flux`, `tabulated_flux` (pinned via `.github/et-legacy-testdata-ref`; only `tabulated_flux` fetches LS220); `PLM_reconstruction`, `WENOZ_reconstruction` |
| `neutrinos` | Table-free physics/classifier fallback plus optically thin gas, constant density sphere, and luminosities |
| `induction-interpolators` / `induction-flux` | Cell/vertex interpolation variants and vector-potential HLL flux variants |
| `compose-regularized-eos` (Ubuntu-GCC only) | Python CompOSE 100% line/branch coverage, synthetic conversion, and StellarCollapse C integration |

### Composite Actions, Coverage Caveats, And `codecov.yml`

- `.github/actions/OS_setup/action.yml` installs compiler and HDF5 dependencies.
- `.github/actions/compile_GRHayL/action.yml` runs `configure`, `make tests datagen`, and `make install` (data generators are compiled only, not executed).
- `.github/actions/code-coverage/action.yml` selects compiler/OS collection commands—Linux GCC uses `gcovr`; Linux clang runs `llvm-profdata` and `llvm-cov` to export `coverage.lcov` (under a comment noting expected coverage files are still not generated); Linux Intel and macOS bodies are commented out—then invokes `codecov/codecov-action@v5` unconditionally whenever called.
- `codecov.yml` defines per-gem components, ignores test paths, targets 100% global patch coverage and 100% zero-threshold CompOSE project/patch coverage, defers notifications until explicit finalization, and documents validation via `curl -X POST --data-binary @codecov.yml https://codecov.io/validate`. The Ubuntu-GCC `compose-regularized-eos` job uploads only `compose-coverage.xml` under flag `compose` with search disabled and fail-on-error enabled.

Claim status: contested; contradiction: CONTR-0006. See [CONTR-0006](../contradictions.md#contr-0006).

Claim status: contested; contradiction: CONTR-0007. See [CONTR-0007](../contradictions.md#contr-0007).

Claim status: contested; contradiction: CONTR-0008. See [CONTR-0008](../contradictions.md#contr-0008).

In `.github/actions/compile_GRHayL/action.yml`, the macOS GCC `CC` command substitution contains an unmatched quote in its `egrep` call, and the macOS clang branch passes `--cflags` twice so `configure` retains only the later include flag and drops coverage flags. In Ubuntu-GCC, `codecov-finalize` uses `always()` and polls recent returned peer runs (macOS GCC, Ubuntu Clang, Ubuntu Intel) for `completed` without requiring every named workflow to appear, succeed, or upload coverage; macOS Clang is excluded.

Evidence (CI behavior): `.github/actions/compile_GRHayL/action.yml::--cflags`; corroboration: `.github/actions/code-coverage/action.yml::codecov/codecov-action@v5` and `github-actions-Ubuntu-gcc.yml::codecov-finalize`.

### Scripted Runner And Read-Only YAML Check

`.github/run_tests.sh` runs `./configure -r`, `make tests datagen`, downloads `GRHayL/TestData` binary reference files and StellarCollapse EOS archives into repo root, executes compiled unit tests (excluding `unit_test_WENOZ_reconstruction`, `unit_test_tabulated_eos_compose`, and `unit_test_con2prim_debug`), invokes `unit_test_code_error` over keys `0`–`90`, runs `pyghl append` on the SLy4 table before NN-enabled tabulated Con2Prim, and cleans tracked downloads on `EXIT`. For full runner behavior and preexisting-file caveats, see [Runner And Build Outputs](../tests/runner-and-build-outputs.md).

Inspect workflow YAML without executing shell commands using PyYAML's `BaseLoader` (which preserves `on` as text instead of YAML 1.1 boolean):

```bash
python3 -c 'import pathlib, yaml
for p in sorted(pathlib.Path(".github/workflows").glob("*.yml")):
    d = yaml.load(p.read_text(), Loader=yaml.BaseLoader)
    print(p, sorted(d["on"]), sorted(d["jobs"]))'
```

Evidence (CI behavior): `run_tests.sh::unit_test_code_error`; corroboration: `et-legacy-testdata-ref`.

## Sources

- [.github/workflows/github-actions-Ubuntu-gcc.yml](../../.github/workflows/github-actions-Ubuntu-gcc.yml) - Ubuntu GCC matrix, `paths-ignore`, `compose-regularized-eos`, `codecov-finalize`
- [.github/workflows/github-actions-Ubuntu-clang.yml](../../.github/workflows/github-actions-Ubuntu-clang.yml) - Ubuntu Clang matrix and `--disable-hdf5` `c2p_nn_guess` leg
- [.github/workflows/github-actions-Ubuntu-intel.yml](../../.github/workflows/github-actions-Ubuntu-intel.yml) - Ubuntu Intel/ICX matrix
- [.github/workflows/github-actions-MacOS-gcc.yml](../../.github/workflows/github-actions-MacOS-gcc.yml) - macOS GCC matrix
- [.github/workflows/github-actions-MacOS-clang.yml](../../.github/workflows/github-actions-MacOS-clang.yml) - macOS Clang matrix
- [.github/actions/OS_setup/action.yml](../../.github/actions/OS_setup/action.yml) - OS dependency installation
- [.github/actions/compile_GRHayL/action.yml](../../.github/actions/compile_GRHayL/action.yml) - configure/compile/install action and macOS quoting/flag issues
- [.github/actions/code-coverage/action.yml](../../.github/actions/code-coverage/action.yml) - compiler-specific coverage collection and Codecov upload
- [.github/run_tests.sh](../../.github/run_tests.sh) - scripted local test replay driver
- [.github/et-legacy-testdata-ref](../../.github/et-legacy-testdata-ref) - pinned TestData reference for flux/ET_Legacy downloads
- [codecov.yml](../../codecov.yml) - coverage components, CompOSE 100% policy, deferred notifications
- [GitHub Actions workflow syntax](https://docs.github.com/en/actions/reference/workflows-and-actions/workflow-syntax#onpushpull_requestpull_request_targetpathspaths-ignore) - `paths-ignore` push/PR semantics

## See Also

- Depends on: [Build And Configure](build-and-configure.md)
- Validated by: [Runner And Build Outputs](../tests/runner-and-build-outputs.md)
- See also: [GRHayLib Verification](../implementations/grhaylib-verification.md)
- See also: [Contradictions](../contradictions.md#contr-0005)
