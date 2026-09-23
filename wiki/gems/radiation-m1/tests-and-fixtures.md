# Radiation M1 Tests And Fixtures

This page routes the grey neutrino M1 test binaries, stored reference data,
generated provider table, CI commands, and evidence limits. The source files,
`Unit_Tests/run_m1_tests.sh`, fixture metadata, `configure`, and CI actions are
the authorities; this page does not replace the numerical contracts or
function-level header documentation.

Read this page with the [Radiation M1 hub](../radiation-m1.md), the
[neutrino M1 contract](neutrino-m1-contract.md), the
[API/build boundary](api-build-boundary.md), and the
[M1 unit-test guide](../../../Unit_Tests/README.m1.md).

## Scoped runner

`Unit_Tests/run_m1_tests.sh` is the canonical scoped runner. It uses the
configured checkout's `build/lib` rather than an installed library and runs
the following ten targets in order:

```text
unit_test_m1_closure_fallback
unit_test_m1_diffusion_flux
unit_test_m1_error_handling
unit_test_m1_fd_jacobian
unit_test_m1_neutrino_rusanov_flux
unit_test_m1_neutrino_seeded_invariants
unit_test_m1_neutrino_source_update
unit_test_m1_rate_provider
unit_test_m1_thcm1_blended_rusanov
unit_test_rusanov_flux
```

After configuring the desired HDF5 mode, use:

```sh
bash Unit_Tests/run_m1_tests.sh --build
```

`--build-only` builds these targets and stops. The runner passes the fixture
directory to the stored-fixture consumers and passes `--generated-fixture` to
the provider test. It does not regenerate the retained reference data.

The dedicated
[Radiation M1 action](../../../.github/actions/run_m1/action.yml) configures a
compiler/HDF5 combination and runs the same `--build` route. The normal
`.github/run_tests.sh` path does not select this scoped runner. The dedicated
action is a test runner only and does not define a separate M1 coverage
threshold, report, artifact, or Codecov policy. Workflow presence is
configured-execution evidence; an observed run must still be reported
separately.

### Evidence vocabulary

Unless a page explicitly identifies an external campaign receipt, an
“agreement” in this inventory means repository-local replay agreement: current
GRHayL is evaluated on the retained inputs and compared with retained THC_M1
endpoints and, where paired inputs exist, their response. This is distinct from
external admission of a fixture package and from host-level or continuum
validation. A local replay agreement must not be read as a new external
admission or as evidence for a complete Cactus/Einstein Toolkit evolution.

## Claim dispositions

The ten targets in the scoped runner have two evidence dispositions. Four
owners combine their existing local checks with the historical package replay:

| Owner | Stored family | Stored claim |
| --- | --- | --- |
| `unit_test_m1_neutrino_seeded_invariants` | Pointwise closure/moments/stress/geometry/speeds, covariant stress-energy, and instantaneous sources | The documented pointwise, CL-04 stress-energy, and instantaneous frozen-rate operations, with named local-policy exceptions kept as separate checks. |
| `unit_test_m1_neutrino_rusanov_flux` | Neutrino and number-current Rusanov | The Rusanov arithmetic for the caller-prepared operands in the fixture. |
| `unit_test_m1_thcm1_blended_rusanov` | Constant-volume and prepared variable-volume transport | The corresponding discrete four-point transport calls and their input/volume conventions. |
| `unit_test_rusanov_flux` | Generic Rusanov | The generic Rusanov call for its serialized E/F operands. |

The remaining six owners are local analytic, invariant, transactional, solver,
or provider checks: `unit_test_m1_closure_fallback`,
`unit_test_m1_diffusion_flux`, `unit_test_m1_error_handling`,
`unit_test_m1_fd_jacobian`, `unit_test_m1_neutrino_source_update`, and
`unit_test_m1_rate_provider`. The diffusion owner additionally replays the
standalone `Jthick` fixture described below; its optional Fick flux helper
remains local-only.
The provider's generated EOS table is test input. The scoped make manifest is
an execution inventory.

## Standalone test coverage

| Test source | Direct behavior exercised | Inputs and fixture route |
| --- | --- | --- |
| [`unit_test_m1_closure_fallback.c`](../../../Unit_Tests/unit_test_m1_closure_fallback.c) | Eulerian Minerbo admissibility fallback: exact symmetry of every published pressure tensor on both the PSD and exact-zero-flux routes, separate `admissibility_fallback_psd` and `admissibility_fallback_zero_flux` accounting, and the documented PSD regime boundary in the aligned, antialigned, and transverse directions. | No stored fixture. The test sweeps flux directions and fluid speeds against a flat metric. |
| [`unit_test_m1_diffusion_flux.c`](../../../Unit_Tests/unit_test_m1_diffusion_flux.c) | `Jthick`, harmonic diffusion coefficient, face-normal proper distance, shared diffusion flux, and the public neutrino diffusion wrapper. It covers active/inactive gates, moving/curved metrics, invalid inputs, overflow rejection, and transactional/no-op behavior. | [`Unit_Tests/data/jthick_thcm1.fixture`](../../../Unit_Tests/data/jthick_thcm1.fixture) stores the six paired `Jthick` inputs/outputs. The optional Fick helper remains local coverage. |
| [`unit_test_m1_error_handling.c`](../../../Unit_Tests/unit_test_m1_error_handling.c) | Shared M1 initialization, closure, repair, stress, diagnostics, geometry/matter sources, fixture parsing, invalid-state/error-code mappings, overflow handling, and unchanged-output rejection. | No external fixture. The test owns deterministic boundary cases. |
| [`unit_test_m1_fd_jacobian.c`](../../../Unit_Tests/unit_test_m1_fd_jacobian.c) | Finite-difference residual/Jacobian construction, public implicit convergence, Newton input and callback validation, admissible projection, retry/backtracking, and transactional failure paths. | No stored fixture. The test supplies frozen rates and local states. |
| [`unit_test_m1_neutrino_rusanov_flux.c`](../../../Unit_Tests/unit_test_m1_neutrino_rusanov_flux.c) | Five-component neutrino Rusanov flux, number-current construction, physical number flux, closure reuse, validation boundaries, seeded paired cases, and stored neutrino/current Rusanov comparisons. | `rusanov_neutrino.dat` and `rusanov_neutrino_current.dat` from `Unit_Tests/data/m1_thcm1/`; `--fixture-dir PATH` selects another fixture directory. |
| [`unit_test_m1_neutrino_seeded_invariants.c`](../../../Unit_Tests/unit_test_m1_neutrino_seeded_invariants.c) | Pointwise closure, comoving moments, stress-energy, geometry sources, wave speeds, source exchange, species and transactional invariants, plus the retained instantaneous source corpus and covariant stress-energy fixtures. | `pointwise_closure_moments.m1`, `stress_energy.m1`, and `m1_thcm1_instantaneous_sources.m1`. The pointwise corpus has 56 paired records; the stress-energy corpus has 1,024 stored pairs; the source corpus has 168 records, including 162 ordinary records and six named local-policy checks. |
| [`unit_test_m1_neutrino_source_update.c`](../../../Unit_Tests/unit_test_m1_neutrino_source_update.c) | Deterministic property cases for explicit thin updates, frozen-rate source updates, diagnostics, mean-energy policy, pair-source updates, terminal/no-update publication, lepton exchange, and invalid inputs. | No stored fixture. The test creates admissible states and rate bundles locally. |
| [`unit_test_m1_rate_provider.c`](../../../Unit_Tests/unit_test_m1_rate_provider.c) | Table-free reference provider, cache reuse and generation changes, channel masks, failure/recovery policies, bounds policies, pair channels, diagnostics, raw-kernel boundaries, and production NRPyLeakage provider behavior when HDF5 is enabled. | `--generated-fixture` creates and removes a deterministic 3 x 3 x 3 StellarCollapse-compatible HDF5 table. An external table path is also accepted by the executable. In no-HDF5 mode the table-backed block is unavailable and the table-free checks remain. |
| [`unit_test_m1_thcm1_blended_rusanov.c`](../../../Unit_Tests/unit_test_m1_thcm1_blended_rusanov.c) | Canonical pointwise four-point blended transport, variable-volume prepared transport, limiter/sawtooth and opacity branches, metric densitization, policy rejection, and transactional validation. | Constant-volume shards `transport_four_point_d0.dat`, `transport_four_point_d1.dat`, and `transport_four_point_d2.dat`; variable-volume shards `transport_four_point_varying_d0.dat`, `transport_four_point_varying_d1.dat`, `transport_four_point_varying_d2.dat`, and `transport_four_point_varying_controls.dat`. The retained corpora contain 57,384 constant-volume pairs and 61,992 variable-volume pairs. |
| [`unit_test_rusanov_flux.c`](../../../Unit_Tests/unit_test_rusanov_flux.c) | Shared component-wise Rusanov arithmetic, M1 physical/Rusanov boundary validation, typed M1 flux cases, and stored generic Rusanov comparisons. | `rusanov_generic.dat`; `--fixture-dir PATH` selects the fixture directory. |

The fixture consumers validate the versioned grammar, record counts, policies,
finiteness, IDs, and the operation-specific checks defined by each adapter.
The stored corpora are consumed by the repository-local tests; the runner does
not regenerate them or download external data.

## Fixture inventory

The tracked fixture directory is
[`Unit_Tests/data/m1_thcm1/`](../../../Unit_Tests/data/m1_thcm1/):

- [`README.md`](../../../Unit_Tests/data/m1_thcm1/README.md) describes the
  pointwise, instantaneous-source, transport, and Rusanov families and their
  provenance boundaries.
- [`TRANSPORT.md`](../../../Unit_Tests/data/m1_thcm1/TRANSPORT.md) documents
  the constant-volume and prepared variable-volume transport corpora,
  direction shards, volume fields, and comparison policy.
- [`SOURCE.md`](../../../Unit_Tests/data/m1_thcm1/SOURCE.md) documents the
  168 instantaneous frozen-rate source records and the six explicit local
  policy classifications.
- [`RUSANOV.md`](../../../Unit_Tests/data/m1_thcm1/RUSANOV.md) documents the
  neutrino, number-current, and generic Rusanov corpora and their caller-input
  boundary.
- [`STRESS_ENERGY.md`](../../../Unit_Tests/data/m1_thcm1/STRESS_ENERGY.md)
  documents the covariant stress-energy corpus, input layout, tensor lowering,
  and fixture boundary.
The runner validates package membership and payload integrity with
`Unit_Tests/data/m1_thcm1/audit_package.py` before replay.

The standalone [`jthick_thcm1.fixture`](../../../Unit_Tests/data/jthick_thcm1.fixture)
is intentionally outside that historical package manifest. It retains the
same generic paired grammar, but only the six existing seeded radiation RNG
records needed by `ghl_m1_compute_Jthick`; it is loaded directly by the
diffusion test and does not add a new package owner or package-level fixture
family.

Offline exporters and retained fixture metadata are regeneration/provenance
tools, not runtime dependencies. Current GRHayL outputs do not replace tracked
fixture values.

## Configuration and evidence limits

- `configure` discovers these as ordinary `unit_test_*.c` targets. The exact
  HDF5 filtering and `GHL_DISABLE_HDF5` behavior comes from `configure`; do
  not infer mode availability from source presence alone.
- The HDF5-enabled provider test reaches the production table-backed provider
  through its generated table. The no-HDF5 route deliberately returns the
  disabled-HDF5 error before the table-backed assembly path.
- The scoped runner audits the stored package before executing its configured
  targets. This checkout contains no multi-toolchain result inventory or
  branch-coverage report; workflow selection is execution configuration, not
  evidence that a remote run has completed.
- These tests establish library-level unit, invariant, transactional, and
  retained-reference behavior. They do not establish Cactus/Einstein Toolkit
  host integration, mesh evolution, reconstruction, AMR, schedules, coupled
  matter stepping, or external physical validation.

## Ground truth references

- [`GRHayL/include/ghl_m1.h`](../../../GRHayL/include/ghl_m1.h)
- [`GRHayL/include/ghl_neutrino_rate_provider.h`](../../../GRHayL/include/ghl_neutrino_rate_provider.h)
- [`GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md`](../../../GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md)
- [`GRHayL/Radiation/TRACEABILITY.md`](../../../GRHayL/Radiation/TRACEABILITY.md)
- [`docs/raw/Radiation.dox`](../../../docs/raw/Radiation.dox)
- [`Unit_Tests/README.m1.md`](../../../Unit_Tests/README.m1.md)
- [`Unit_Tests/run_m1_tests.sh`](../../../Unit_Tests/run_m1_tests.sh)
- [`wiki/test-map.md`](../../test-map.md)
