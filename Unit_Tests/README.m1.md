# Radiation M1 unit tests

The two-layer evidence boundary for this directory is summarized in
[`README.md`](README.md). This page records the current scoped M1 inventory,
stored fixture families, and their evidence limits.

The existing M1 test executables exercise shared Radiation kernels and grey
three-species neutrino operations. Their randomized, analytic, invariant, and
transactional checks remain part of normal execution. Stored THC_M1 comparisons
are incorporated into the relevant existing executables; the historical package
provenance and operation-specific coverage are documented with the data in
`data/m1_thcm1/`. The standalone `data/jthick_thcm1.fixture` is the direct
paired replay for the thick-limit scalar. Its six baseline/perturbed pairs all
agree with the THC_M1 thick-pressure comoving-energy projection, including all
six paired responses; the optional Fick flux remains local-only.

Local regressions in `unit_test_m1_neutrino_seeded_invariants` cover nonzero
fluxes whose squared norms underflow or whose `E/F^2` intermediates overflow,
while retaining the distinct exact-zero closure policy. The
`unit_test_m1_neutrino_source_update` regressions check mean-energy bounds at
the final paired endpoint, including valid endpoints reached through an
out-of-bounds intermediate state and transactional rejection of invalid final
endpoints.

The closure regression checks that normalized pressure and closure scalars
remain unchanged when energy scaling crosses the `E^2` overflow boundary in
a nonidentity metric. Source-update regressions check charged-current lepton
conservation at stiff backward-Euler endpoints, separately from number-floor
repairs and the optional equilibrium-number projection.

## Build and run

From the repository root, for a build without HDF5:

```sh
./configure --noomp --disable-hdf5 --prefix="$PWD"
bash Unit_Tests/run_m1_tests.sh --build
```

For HDF5 support, omit `--disable-hdf5` and configure its installation normally.
Use a fresh build when changing compiler or configuration. The runner builds
only the named existing M1 targets and their required library/test-support
dependencies. It does not build all test executables or an Einstein Toolkit
configuration. `--build-only` builds the targets without running them. With no
option, the runner executes already-built tests. It resolves the repository root
from its location, so it can also be invoked from another working directory.
On a clean checkout it expands the tracked `data/m1_thcm1/payloads.tar.gz`
before auditing the retained fixtures. Existing plaintext payloads are audited
without replacement; no external download is needed.

The runner uses `build/lib` for the current checkout's shared library. It passes
`--generated-fixture` to the rate-provider test: HDF5 builds create and remove a
deterministic test-local EOS table and exercise the table-backed provider.
Without HDF5, the same invocation runs the available table-free checks. The
provider executable retains its optional external EOS-table argument. This
table is provider-test input, not an independent reference output.

## NRPyLeakage adapter scope

The M1 NRPyLeakage backend uses the NRPyLeakage-owned, source-private helper
headers for local nucleon blocking, reaction-shifted beta moments, and
bremsstrahlung moments:

- [`NRPyLeakage_nucleon_blocking.h`](../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_nucleon_blocking.h)
- [`NRPyLeakage_rate_helpers.h`](../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_rate_helpers.h)

Direct source-private cross-gem inclusion is an approved boundary for this
adapter. The M1 manifest records the dependency; the helper implementations
remain owned by NRPyLeakage, are not copied into Radiation, and are not
installed as public API. M1 remains tau-free: optical-depth suppression,
leakage luminosity/source assembly, and legacy leakage finite-output fallback
policies do not enter the M1 rate bundle. M1 keeps its own validation and
provider recovery behavior.

## Stored reference boundary

Normal fixture replay requires neither THC_M1 source nor `THCM1_ROOT` nor the
external Verification directory or an external reference-data download.
Reference generation is an offline operation. The stored package is not a live
campaign runner.
Stored expected values must come from a retained THC producer record with
reviewed provenance and must remain paired with the actual physical inputs it
consumed. Current GRHayL
outputs are computed during testing and never used to regenerate expectations.

The fixture-family documents record the retained producer receipt, sidecar,
command capture, consumed-input stream, and selected result evidence required
by each offline exporter. They document provenance for the repository-local
package; `wiki/gems/radiation-m1/tests-and-fixtures.md` is the status authority
for this checkout. External `Verification/` admission manifests and mutable
campaign results are neither required nor authoritative; when used as a
source/evaluator build input, they remain provenance only.

The exact variable-volume transport family is loaded by the existing
`unit_test_m1_thcm1_blended_rusanov` executable from
`transport_four_point_varying_d0.dat`, `transport_four_point_varying_d1.dat`,
`transport_four_point_varying_d2.dat`, and the
`transport_four_point_varying_controls.dat` control shard. It uses the existing 50-field
transport input layout: face volume, speed, opacity, spacing, theta, minimum
dissipation, four cell volumes, four five-component stencil states, and two
five-component physical face fluxes. The test adapter validates the volumes
and prepares `V_cell U` and `V_face F` before calling the required
`ghl_m1_compute_neutrino_four_point_volume_weighted_transport_flux`
transport API. Missing variable-volume shards are a hard test failure; no
retained reference output is generated or substituted at runtime.

Reference comparisons supplement local tests. Matching a pointwise or prepared
face output does not establish full-grid evolution, continuum convergence, or
whole-code equivalence. Operations with different numerical or physical models
retain local tests and an explicit comparison limitation. The discrete-operator
equivalence workflow, live THC_M1 execution, host evolution, and provider/rate
workflows remain separate claims.

Five of the nine executables replay stored baseline/perturbed THC
outputs: `unit_test_m1_neutrino_seeded_invariants`,
`unit_test_m1_neutrino_rusanov_flux`, `unit_test_m1_thcm1_blended_rusanov`, and
`unit_test_rusanov_flux` use the historical package, while
`unit_test_m1_diffusion_flux` uses the standalone six-pair `Jthick` fixture.
The established corpora have these distinctions:

| Stored family | Pairs | Changed-input pairs | Identical-input controls |
| --- | ---: | ---: | ---: |
| Pointwise | 56 | 31 | 25 |
| Instantaneous sources | 168 | 120 | 48 |
| Constant-volume transport | 57,384 | 38,880 | 18,504 |
| Varying-volume transport, including control shard | 61,992 | 38,880 | 23,112 |
| Neutrino Rusanov, original corpus | 1,024 | 1,024 | 0 |
| Neutrino Rusanov, supplemental nonzero current | 1,024 | 1,024 | 0 |
| Generic Rusanov | 1,024 | 683 | 341 |
| Stress-energy tensor | 1,024 | 1,024 | 0 |
| Thick-limit `Jthick` | 6 | 6 | 0 |

Pointwise and instantaneous-source families include respectively two and six
named exact-zero-flux moving-fluid policy differences. Those records exercise
local policy checks; they are not THC agreement claims. Identical-input records
are controls, and changed inputs need not change every output. Rusanov replay
evaluates current GRHayL at the baseline input; its retained THC perturbation
supplies the response envelope and is not recomputed by GRHayL. Pointwise,
instantaneous-source, stress-energy, and prepared transport are paired
replays: current GRHayL is evaluated at both paired inputs, and the current
baseline/perturbed response is compared with the retained THC response. The
prepared face operands still limit Rusanov/transport agreement to the discrete
operation, rather than independently validating physical flux and closure
preparation.

For the pointwise owner, ordinary records evaluate current GRHayL at both
baseline and perturbed inputs, and compare the current response with the paired
THC_M1 response. The two named zero-flux policy records remain explicit local
two-state admissibility checks.

The seeded-invariants owner also replays `stress_energy.m1`. Its complete
21-value state/ADM/velocity input is evaluated through the public neutrino
stress-energy wrapper at both paired endpoints, lowered with the full ADM
four-metric, and compared with the retained THC_M1 `assemble_rT` endpoints and
response. This remains a discrete stress-energy check; it does not claim host
evolution or continuum agreement.

There is no stored THC finite-step pair/implicit endpoint corpus. Independent
local collision-model oracles and solver checks cover those operations. The
explicit source-update, repair, and failure-path checks therefore remain local
evidence. See
[RUSANOV.md](data/m1_thcm1/RUSANOV.md) for the supplemental nonzero-current
fixture family, and [STRESS_ENERGY.md](data/m1_thcm1/STRESS_ENERGY.md) for the
covariant stress-energy corpus. The runner validates the retained package with
`audit_package.py` before executing its consumers.

## CI

The broad `.github/run_tests.sh` runner does not invoke the scoped M1 suite.
The compiler/OS workflows select dedicated Radiation jobs through
`.github/actions/run_m1`; these configure the selected compiler and HDF5 mode,
then build and execute the same scoped M1 targets. The Ubuntu GCC Radiation
jobs compile with gcov flags and invoke the shared GRHayL coverage action after
the scoped run in both HDF5 modes. That action uploads a gcovr Cobertura report
filtered to `GRHayL/Radiation/`. The `radiation_m1` Codecov component has a
100% project coverage target for that path. The M1 action defines no
changed-path gate. Workflow selection is not a claim that a remote CI run has
already passed.
