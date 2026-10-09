# Radiation M1 Tests And Fixtures

This page routes the grey neutrino M1 test binaries, stored reference data,
generated provider table, CI commands, and evidence limits. The source files,
`.github/run_tests.sh`, fixture metadata, `configure`, and CI actions are
the authorities; this page does not replace the numerical contracts or
function-level header documentation.

Read this page with the [Radiation M1 hub](../radiation-m1.md), the
[neutrino M1 contract](neutrino-m1-contract.md), the
[API/build boundary](api-build-boundary.md), and the
[M1 unit-test guide](../../../docs/raw/Radiation_unit_tests.md). The
[`run_m1` composite action](../../../.github/actions/run_m1/action.yml) replays the
pinned TestData publication by default: it retries the ordinary
[`.github/radiation-testdata-ref`](../../../.github/radiation-testdata-ref)
acquisition and forwards its optional `fixture-dir` input to the runner as
`M1_FIXTURE_DIR`. A nonempty
`fixture-dir` input instead replays a supplied raw or gzip directory and skips
downloads. The ordinary runner with unset `M1_FIXTURE_DIR` downloads and
prepares every pinned member by default; individual executable invocations
without fixture arguments remain local-only.

## Build and execution

`configure` discovers the ordinary `unit_test_*.c` targets. Build only the
affected `test/unit_test_*` targets and their required library dependencies;
the [unit-test guide](../../../docs/raw/Radiation_unit_tests.md#build-and-run)
lists the M1 scope. No Einstein Toolkit build is needed.

The single [runner](../../../.github/run_tests.sh) takes
`[all|m1] [configure arguments...]`: `all` (the default) executes the M1 tests
alongside the existing suite, and `m1` builds with `make tests` and runs only
the Radiation M1 tests. Each compiler workflow in
[`.github/workflows/`](../../../.github/workflows/) has a `radiation-m1` job
that calls the
[`run_m1` action](../../../.github/actions/run_m1/action.yml), which selects the
compiler and HDF5 mode and ends with `run_tests.sh m1`, in independent
HDF5-enabled and disabled builds. Only the
[Ubuntu GCC workflow](../../../.github/workflows/github-actions-Ubuntu-gcc.yml)
job builds with coverage flags and then runs the
[M1 coverage action](../../../.github/actions/m1-code-coverage/action.yml), which
generates and uploads the gcovr report for `GRHayL/Radiation/`, including
executable private headers, and enforces its gates; exemptions require
independently checked proofs. The other `radiation-m1` jobs are plain test runs.
See the [unit-test guide](../../../docs/raw/Radiation_unit_tests.md#ci) for
measurement and provider evidence boundaries. Local checks need no stored
fixture publication. Workflow configuration and local runner verification are
separate from an observed remote CI pass.

Stored-reference consumers accept `--fixture-dir PATH`; the diffusion owner
accepts `--fixture PATH` for Jthick. The provider's `--generated-fixture` creates
local EOS input and does not generate trusted expected M1 outputs.

The TestData storage form is individual `.bin.gz` members. The
[preparation helper](../../../.github/prepare_m1_fixtures.sh) restores the
original `.bin` layout into an empty caller-owned directory before replay;
raw members can also be copied through that interface. The ordinary runner
with unset `M1_FIXTURE_DIR` downloads every named `radiation/*.bin.gz` file
from raw.githubusercontent.com at the full commit in
`.github/radiation-testdata-ref` through its existing shared download
machinery, prepares the raw members in a private `${TMPDIR:-/tmp}` directory
registered with its EXIT cleanup, and replays them by default. Setting
`M1_FIXTURE_DIR` to a supplied raw or gzip package bypasses downloads, leaves
caller inputs unchanged, and retains replay: a complete raw directory is
replayed directly, while compressed members are prepared in a temporary
directory. For a sibling local TestData checkout, set
`M1_FIXTURE_DIR="$(pwd)/../TestData/radiation"` when invoking the existing
`run_tests.sh m1` runner from the GRHayL root. This local path requires no
publication or external revision pin.
The [standalone instructions](../../../docs/raw/Radiation_unit_tests.md#build-and-run)
show preparation, executable arguments, and temporary-directory cleanup.

### Evidence vocabulary

The [fresh candidate production record](../../../docs/raw/Radiation_fresh_fixture_production.md)
describes new THC_M1 runs, receipts, consumed inputs, raw outputs, source
snapshots, and exporters for every payload. The accepted package is published
as individual `radiation/*.bin.gz` members in `GRHayL/TestData` and is
downloaded and replayed by default; the
[TestData provenance](https://github.com/GRHayL/TestData/blob/main/radiation/PROVENANCE.md)
describes the published fixture families and binary layout. The imported
[TestData revision file](../../../.github/radiation-testdata-ref)
is the pinned acquisition source for every default runner and action replay.
An explicitly
supplied older package remains historical replay only. Neither fixture replay
nor the fresh producer run proves host-level evolution or continuum validation.

## Claim dispositions

These scoped owners combine their existing local checks with the historical
package replay:

| Owner | Stored family | Stored claim |
| --- | --- | --- |
| `unit_test_m1_neutrino_seeded_invariants` | Pointwise closure/moments/stress/geometry/speeds, covariant stress-energy, and instantaneous sources | The documented pointwise, CL-04 stress-energy, and instantaneous frozen-rate operations, with named local-policy exceptions kept as separate checks. |
| `unit_test_m1_neutrino_rusanov_flux` | Neutrino and number-current Rusanov | The Rusanov arithmetic for the caller-prepared operands in the fixture. |
| `unit_test_m1_thcm1_blended_rusanov` | Constant-volume and prepared variable-volume transport | The corresponding discrete four-point transport calls and their input/volume conventions. |
| `unit_test_rusanov_flux` | Generic Rusanov | The generic Rusanov call for its serialized E/F operands. |

The other scoped owners run local analytic, invariant, transactional, solver,
or provider checks: `unit_test_m1_closure_fallback`,
`unit_test_m1_diffusion_flux`, `unit_test_m1_error_handling`,
`unit_test_m1_fd_jacobian`, `unit_test_m1_neutrino_source_update`, and
`unit_test_m1_rate_provider`. The diffusion owner additionally replays the
standalone `Jthick` fixture described below; its optional Fick flux helper
remains local-only.
The provider's generated EOS table is test input. The configured Makefile registers the discovered test targets.

## Standalone test coverage

| Test source | Direct behavior exercised | Inputs and fixture route |
| --- | --- | --- |
| [`unit_test_m1_closure_fallback.c`](../../../Unit_Tests/unit_test_m1_closure_fallback.c) | Eulerian Minerbo admissibility fallback: exact symmetry of every published pressure tensor on the fallback as well as the primary path, publication through both the PSD and exact-zero-flux fallback causes, and the documented PSD regime boundary in the aligned, antialigned, and transverse directions. | No stored fixture. The test sweeps flux directions and fluid speeds against a flat metric. |
| [`unit_test_m1_diffusion_flux.c`](../../../Unit_Tests/unit_test_m1_diffusion_flux.c) | `Jthick`, harmonic diffusion coefficient, face-normal proper distance, shared diffusion flux, and the public neutrino diffusion wrapper. It covers active/inactive gates, moving/curved metrics, invalid inputs, overflow rejection, and transactional/no-op behavior. | External `jthick_thcm1.bin` stores the paired `Jthick` inputs/outputs; `--fixture PATH` selects it and the ordinary runner supplies it under `M1_FIXTURE_DIR`. No text copy is tracked in this checkout. The optional Fick helper remains local coverage. |
| [`unit_test_m1_error_handling.c`](../../../Unit_Tests/unit_test_m1_error_handling.c) | Shared M1 initialization, closure, repair, stress, diagnostics, geometry/matter sources, fixture parsing, invalid-state/error-code mappings, overflow handling, and unchanged-output rejection. | No external fixture. The test owns deterministic boundary cases. |
| [`unit_test_m1_fd_jacobian.c`](../../../Unit_Tests/unit_test_m1_fd_jacobian.c) | Finite-difference residual/Jacobian construction, public implicit convergence, Newton input and callback validation, admissible projection, retry/backtracking, small-correction no-root and energy-floor rejection, and transactional failure paths. | No stored fixture. The test supplies frozen rates and local states. |
| [`unit_test_m1_neutrino_rusanov_flux.c`](../../../Unit_Tests/unit_test_m1_neutrino_rusanov_flux.c) | Five-component neutrino Rusanov flux, number-current construction, physical number flux, closure reuse, validation boundaries, seeded paired cases, and stored neutrino/current Rusanov comparisons. | External retained `rusanov_neutrino.bin` and `rusanov_neutrino_current.bin`; `--fixture-dir PATH` selects their directory. |
| [`unit_test_m1_neutrino_seeded_invariants.c`](../../../Unit_Tests/unit_test_m1_neutrino_seeded_invariants.c) | Pointwise closure, comoving moments, stress-energy, geometry sources, wave speeds, source exchange, species and transactional invariants, plus the retained instantaneous source corpus and covariant stress-energy fixtures. This owner holds the exact-zero and seeded-state cases; binary64 provider arithmetic underflow/overflow rejection belongs to the rate-provider owner below. Local checks cover `ghl_m1_compute_neutrino_explicit_rhs_sources` geometry-only and interaction-enabled composition, optional number output, and unchanged outputs on failure. | `pointwise_closure_moments.bin`, `stress_energy.bin`, and `m1_thcm1_instantaneous_sources.bin`. The retained pointwise, stress-energy, and source corpus formats and record extents are defined by the adapter schemas and the [fixture-family documents](../../../docs/raw/m1_thcm1/README.md), including the named local-policy checks. The explicit RHS checks are local, not retained THC_M1 endpoint comparisons. |
| [`unit_test_m1_neutrino_source_update.c`](../../../Unit_Tests/unit_test_m1_neutrino_source_update.c) | Deterministic property cases for explicit thin updates, frozen-rate source updates, diagnostics, mean-energy policy, pair-source updates, terminal/no-update publication, lepton exchange, and invalid inputs. | No stored fixture. The test creates admissible states and rate bundles locally. |
| [`unit_test_m1_rate_provider.c`](../../../Unit_Tests/unit_test_m1_rate_provider.c) | Test-local synthetic reference checks plus independent production-provider cache, channel, bounds, physical recovery, diagnostics and raw-kernel checks, including binary64 provider arithmetic underflow/overflow rejection. | `--generated-fixture` creates and removes a deterministic StellarCollapse-compatible HDF5 table whose grid extents the executable defines. An external table path is also accepted by the executable. In no-HDF5 mode the table-backed block is unavailable and test-local reference checks and production disabled-HDF5 transactions remain. |
| [`unit_test_m1_thcm1_blended_rusanov.c`](../../../Unit_Tests/unit_test_m1_thcm1_blended_rusanov.c) | Canonical pointwise four-point blended transport, variable-volume prepared transport, limiter/sawtooth and opacity branches, metric densitization, policy rejection, and transactional validation. | Constant-volume shards `transport_four_point_d0.bin`, `transport_four_point_d1.bin`, and `transport_four_point_d2.bin`; variable-volume shards `transport_four_point_varying_d0.bin`, `transport_four_point_varying_d1.bin`, `transport_four_point_varying_d2.bin`, and `transport_four_point_varying_controls.bin`. The retained corpus payload sizes are defined by the
[transport fixture document](../../../docs/raw/m1_thcm1/TRANSPORT.md). |
| [`unit_test_rusanov_flux.c`](../../../Unit_Tests/unit_test_rusanov_flux.c) | Shared component-wise Rusanov arithmetic, M1 physical/Rusanov boundary validation, typed M1 flux cases, and stored generic Rusanov comparisons. | `rusanov_generic.bin`; `--fixture-dir PATH` selects the fixture directory. |

The fixture consumers validate the versioned binary layout, record counts, policies,
finiteness, IDs, and the operation-specific checks defined by each adapter.
The shared reader also accepts the original text layout for the repository-local
Jthick fixture and parser rejection checks.
Case and pair ID uniqueness is checked through sorted borrowed IDs without
changing replay order. The Rusanov consumers evaluate both input roles, retain
the baseline-envelope gate, and additionally enforce the named strict paired
endpoint/response policy. Its propagated response bound checks consistency of
the endpoints and their difference, not independent derivative accuracy.
Requested replay reads external retained data; it never regenerates expected
values from current GRHayL outputs.

## Fixture inventory

The large retained payload archive is no longer part of this repository. A
replacement candidate with all pointwise, instantaneous-source, stress-energy,
Rusanov, constant-volume transport, variable-volume transport, and Jthick
families has been generated from fresh THC_M1 executions; its accepted gzip
payloads are published as `radiation/*.bin.gz` members in `GRHayL/TestData`,
and its producer-evidence
archive is retained outside the repository. The published fixture families
and binary layout are described by the
[TestData provenance](https://github.com/GRHayL/TestData/blob/main/radiation/PROVENANCE.md),
and the operation formats and claim boundaries are documented in:

- [Fixture overview](../../../docs/raw/m1_thcm1/README.md)
- [Transport](../../../docs/raw/m1_thcm1/TRANSPORT.md)
- [Instantaneous sources](../../../docs/raw/m1_thcm1/SOURCE.md)
- [Rusanov](../../../docs/raw/m1_thcm1/RUSANOV.md)
- [Stress-energy](../../../docs/raw/m1_thcm1/STRESS_ENERGY.md)

A sibling historical copy is not producer evidence. The fresh candidate has
captured commands, receipts, sidecars, source snapshots, consumed inputs, and
raw THC_M1 outputs for every family. The CI action and the ordinary runner
acquire the pinned TestData publication by default; local invocations can
supply a fixture directory explicitly.
Local replay alone cannot establish
admission. Offline exporters are maintenance tools, never normal test
dependencies.

## Configuration and evidence limits

- `configure` discovers these as ordinary `unit_test_*.c` targets. The exact
  HDF5 filtering and `GHL_DISABLE_HDF5` behavior comes from `configure`; do
  not infer mode availability from source presence alone.
- The HDF5-enabled provider test reaches the production table-backed provider
  through its generated table. The no-HDF5 route deliberately returns the
  disabled-HDF5 error before the table-backed assembly path.
- Explicit replay arguments select retained-reference checks; requested missing
  or incomplete fixtures fail. Workflow selection is execution configuration,
  not evidence that a remote run has completed.
- These tests establish library-level unit, invariant, transactional, and
  retained-reference behavior. They do not establish Cactus/Einstein Toolkit
  host integration, mesh evolution, reconstruction, AMR, schedules, coupled
  matter stepping, or external physical validation.

## Ground Truth

- [`GRHayL/include/ghl_m1.h`](../../../GRHayL/include/ghl_m1.h)
- [`GRHayL/include/ghl_neutrino_rate_provider.h`](../../../GRHayL/include/ghl_neutrino_rate_provider.h)
- [`docs/raw/Radiation_integration_contract.md`](../../../docs/raw/Radiation_integration_contract.md)
- [`docs/raw/Radiation_traceability.md`](../../../docs/raw/Radiation_traceability.md)
- [`docs/raw/Radiation.dox`](../../../docs/raw/Radiation.dox)
- [`docs/raw/Radiation_unit_tests.md`](../../../docs/raw/Radiation_unit_tests.md)
- [`.github/run_tests.sh`](../../../.github/run_tests.sh)
- [`wiki/test-map.md`](../../test-map.md)
