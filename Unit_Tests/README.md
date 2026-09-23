# Unit-test boundary

This page defines the evidence boundary for the focused Radiation M1 unit-test
inventory. The broader test inventory and legacy fixture conventions are routed
through [`wiki/test-map.md`](../wiki/test-map.md) and the repository
[`README.md`](../README.md). The M1-specific details are in
[`README.m1.md`](README.m1.md).

## Two M1 evidence layers

### Local analytic, invariant, and transactional tests

Local tests construct their cases in the executable or in test-local helpers.
They cover identities, limiting cases, admissibility and error paths,
transactional/no-op behavior, solver properties, and provider behavior. A
locally generated state, output, or EOS table is test input or a local oracle;
it is not an independent cross-code reference.

The current scoped runner has these local-only owners:

- `unit_test_m1_diffusion_flux`
- `unit_test_m1_error_handling`
- `unit_test_m1_fd_jacobian`
- `unit_test_m1_neutrino_source_update`
- `unit_test_m1_rate_provider`

The four stored-reference owners also retain local checks around their replay.
Those local checks remain local evidence even when the same executable consumes
a stored fixture.

In particular, `--generated-fixture` on the rate-provider test creates a
deterministic table for provider coverage. It does not generate a trusted M1
output and does not establish a provider or rate cross-code claim.

### Offline stored THC_M1 replays

The repository-local package is
[`data/m1_thcm1/`](data/m1_thcm1/). The current stored-reference owners are:

| Owner | Stored operation families | Claim boundary |
| --- | --- | --- |
| `unit_test_m1_neutrino_seeded_invariants` | Pointwise closure/moments/stress/geometry/speeds and instantaneous frozen-rate sources | The documented pointwise and source operations, including their named local-policy exceptions. |
| `unit_test_m1_neutrino_rusanov_flux` | Neutrino and number-current Rusanov fixtures | The Rusanov operation with the caller-prepared operands recorded by the fixture. |
| `unit_test_m1_thcm1_blended_rusanov` | Constant-volume and variable-volume four-point transport | The corresponding pointwise or prepared-face discrete transport operation, including paired baseline/perturbed response. |
| `unit_test_rusanov_flux` | Generic Rusanov fixture | The generic Rusanov operation with its supplied E/F operands. |

Each stored record retains a baseline/perturbed input pair and the associated
THC_M1 outputs. Pointwise, source, Rusanov, and stress-energy adapters evaluate
current GRHayL at the baseline input, compare that output with the retained
baseline, and use the retained perturbation response as the error envelope.
The prepared-transport adapter evaluates current GRHayL at both inputs and
compares its baseline/perturbed response with the retained THC response. The
named local two-state policy records are the other exceptions. A replay
therefore proves only the named library operation and its serialized input
convention; it does not validate upstream preparation that the fixture stores
as caller input.

Stored expected values are produced offline by the retained THC_M1 producer.
Normal tests read only the tracked package and never invoke THC_M1,
`THCM1_ROOT`, `Verification/`, or an external reference-data download. Current
GRHayL output is never used to regenerate an expected value during a test.

## Fixture provenance and admission

The fixture README and the family documents in `data/m1_thcm1/` describe the
retained producer evidence: completed campaign receipts, canonical sidecars,
recorded commands, consumed-input streams, and the selected result rows where
those artifacts are retained. The offline exporters are maintenance and
regeneration tools; they are not runtime dependencies.

The package is repository-local, but Unit_Tests does not perform live campaign
execution. A passing replay is therefore not a fresh live THC_M1 run; the
current-PR route contracts, audits, and radiation validation page establish the
status for this checkout. External `Verification/` manifests and mutable
campaign results are provenance or tooling inputs only, not an authority for
admission or promotion. Missing operands must remain blocked rather than being
filled with current GRHayL output. The retained provenance also does not claim
that current GRHayL is identical to the historical producer source snapshot.

## Separate workflows and evidence limits

The following remain outside this unit-test claim boundary:

- live THC_M1 execution and mutable `Verification/` campaign results;
- host evolution, mesh loops, reconstruction, AMR, schedules, and downstream
  Einstein Toolkit integration;
- discrete-operator equivalence workflows;
- provider/rate validation as an independent cross-code claim; and
- full-grid, continuum-convergence, physical-validation, or whole-code
  equivalence claims.

No stored finite-step or implicit-endpoint THC_M1 corpus is shipped. The local
source-update, finite-difference/Jacobian, repair, and failure-path tests do not
acquire a stored-reference claim by proximity to a stored fixture. Surfaces not
named in the current stored-owner table remain local or separate until an
admitted fixture and an executable owner exist.
