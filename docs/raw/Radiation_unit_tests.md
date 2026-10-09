# Radiation M1 unit tests

This page defines the evidence boundary for the focused Radiation M1 unit-test
inventory, scoped builds, the stored fixture families, and their limits. The
broader test inventory and legacy fixture conventions are routed through
`wiki/test-map.md` and the repository `README.md`. Fixture-family documents for
the stored THC_M1 package are in the `m1_thcm1/` directory beside this page.

## Two M1 evidence layers

### Local analytic, invariant, and transactional tests

Local tests construct their cases in the executable or in test-local helpers.
They cover identities, limiting cases, admissibility and error paths,
transactional/no-op behavior, solver properties, and provider behavior. A
locally generated state, output, or EOS table is test input or a local oracle;
it is not an independent cross-code reference.

These owners provide local-only checks:

- `unit_test_m1_closure_fallback`
- `unit_test_m1_error_handling`
- `unit_test_m1_fd_jacobian`
- `unit_test_m1_neutrino_source_update`
- `unit_test_m1_rate_provider`

The five stored-reference owners also retain local checks around their replay.
Those local checks remain local evidence even when the same executable consumes
a stored fixture.

In particular, `--generated-fixture` on the rate-provider test creates a
deterministic table for provider coverage. It does not generate a trusted M1
output and does not establish a provider or rate cross-code claim.

### Offline stored THC_M1 replays

Retained data is delivered separately from the library source. Replay uses
`--fixture-dir PATH`, or `--fixture PATH` for Jthick. The stored-reference
owners are:

| Owner | Stored operation families | Claim boundary |
| --- | --- | --- |
| `unit_test_m1_neutrino_seeded_invariants` | Pointwise closure/moments/stress/geometry/speeds, instantaneous frozen-rate sources, and the stress-energy tensor | The documented pointwise, source, and stress-energy operations, including their named local-policy exceptions. |
| `unit_test_m1_neutrino_rusanov_flux` | Neutrino and number-current Rusanov fixtures | The Rusanov operation with the caller-prepared operands recorded by the fixture. |
| `unit_test_m1_thcm1_blended_rusanov` | Constant-volume and variable-volume four-point transport | The corresponding pointwise or prepared-face discrete transport operation, including paired baseline/perturbed response. |
| `unit_test_rusanov_flux` | Generic Rusanov fixture | The generic Rusanov operation with its supplied E/F operands. |
| `unit_test_m1_diffusion_flux` | Standalone thick-limit `Jthick` fixture | The `Jthick` operation at the six stored baseline/perturbed pairs. |

Each stored record retains a baseline/perturbed input pair and THC_M1 outputs.
For the fresh candidate, the
[production record](Radiation_fresh_fixture_production.md)
identifies the producer records and exporters; the removed historical package
had unverified producer attribution. Rusanov adapters evaluate current GRHayL
at both paired inputs. They preserve the baseline comparison using the retained
perturbation response as its error envelope and additionally require strict
endpoint and propagated-response agreement under the fixture's named policy.
The response bound propagates endpoint uncertainty; it is not an independent
derivative-accuracy test. Pointwise, instantaneous-source,
stress-energy, and prepared-transport adapters evaluate current GRHayL at both
paired inputs and compare the current baseline/perturbed response with the
retained THC response. The named local two-state policy records are local
policy checks rather than THC agreement claims. A replay
therefore proves only the named library operation and its serialized input
convention; it does not validate upstream preparation that the fixture stores
as caller input.

The fresh candidate maps all 14 payloads to new THC_M1 producer runs. Replay
reads the requested package and never invokes THC_M1,
`THCM1_ROOT`, or offline producer campaigns to generate expected values. Current
GRHayL output is never used to regenerate an expected value during a test.

## Fixture provenance and admission

The [fixture overview](m1_thcm1/README.md) and family documents describe the
formats and claim boundaries. The
[production record](Radiation_fresh_fixture_production.md) records the
accepted payloads, their replay route, and the maintained
[TestData provenance](https://github.com/GRHayL/TestData/blob/main/radiation/PROVENANCE.md)
that describes the published fixture families and binary layout. The old
payload archive, manifest, and auditor have been removed from the library
source. The published package in
`GRHayL/TestData` holds only the fourteen `radiation/*.bin.gz` members and
`PROVENANCE.md`; the production record and evidence archive retain the family
manifest, receipts, raw producer outputs, source
snapshots, and exporters outside TestData. The ordinary runner downloads and
prepares the fourteen pinned members at the
`.github/radiation-testdata-ref` revision by default. A nonempty
`fixture-dir` (or `M1_FIXTURE_DIR`) instead replays a supplied local raw or
gzip package and skips downloads. A passing
local replay alone is not fresh THC_M1 execution. Missing producer operands
must never be replaced with current GRHayL output.

The shared `ghl_pert_test_fail_with_tolerance` comparison is used for the
retained baseline/response policy after exact scaling when representable. The
pointwise and strict paired fixture policies use an inclusive, input-normalized
bound (the pointwise bound also has a symmetric relative term and a floor),
which that helper cannot express. Their local comparator preserves the stored
policy. The baseline/response comparator uses its own range fallback only when
scaling would overflow or erase a nonzero value.

## Scoped inventory and local regressions

The existing M1 test executables exercise shared Radiation kernels and grey
three-species neutrino operations. Their randomized, analytic, invariant, and
transactional checks remain part of normal execution. Stored THC_M1 comparisons
are incorporated into the relevant existing executables; the historical package
provenance and operation-specific coverage are documented in `docs/raw/m1_thcm1/`.
The standalone Jthick payload is the direct paired replay for the thick-limit
scalar. The existing comparisons match its six baseline/perturbed retained
values and paired responses. The producer-to-payload binding remains
unverified, so this is not an admitted THC_M1 comparison; the optional Fick
flux remains local-only.

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

The generic Rusanov boundary checks retain finite results when an intermediate
dissipation product overflows and when cancellation requires a half-scaled
retry. Configure with `--cflags='-ffp-contract=off'` to verify that correctness
does not depend on the compiler implicitly fusing a multiply and subtraction.
Truly overflowing candidates are rejected without modifying the output buffer.

### Double arithmetic and demonstrated callers

Production Radiation arithmetic uses `double`, not platform-dependent
`long double` or `LDBL_*` acceptance gates. Normalized metric/closure checks
use `DBL_EPSILON`-scaled tolerances and reject inputs outside their accuracy
domain. Range-preserving arithmetic is restricted to the following existing
caller needs; it does not relax realizability, convergence, or endpoint gates.

| Retained operation | Production need | Existing regression owner |
| --- | --- | --- |
| Overflow-safe average, explicit FMA and half-scaled Rusanov retry in `ghl_calculate_Rusanov_flux.c` | Two-state E/F and number wrappers and the canonical four-point core share the scalar expression; prepared face operands carry host volume factors. A representable result must not be rejected solely because its state jump, dissipation product, or rounded flux average overflows. The direct finite path retains its ordinary rounding. | `unit_test_rusanov_flux.c` checks the overflowing-jump result, fused overflowing-product result, half-scaled finite-jump result, and unchanged output on genuine overflow; `unit_test_m1_thcm1_blended_rusanov.c` checks the production four-point subnormal-speed path. |
| Exponent-aligned difference of products in `ghl_m1_utils.h` | Physical E/F and number fluxes subtract `alpha*current - beta*state`; both products can overflow while the coordinate flux is finite. | `check_private_norm_boundaries` in `unit_test_m1_error_handling.c` calls the public physical-flux kernel with overflowing products and checks its finite energy/momentum outputs exactly. |
| Explicit closure FMA and energy-normalized equations | Boosted isotropic radiation cancels large terms in the thick tensor; energy normalization avoids squared-energy overflow without choosing extended precision. | `check_boosted_isotropic_endpoint` in `unit_test_m1_closure_fallback.c` checks the thick endpoint across boosts and energy scales with `DBL_EPSILON`-scaled gates. |
| Double-mantissa/exponent positive products, sums, ratios and comparisons in `ghl_m1_neutrino_implicit.h` | Backward-Euler emission/absorption ratios, charged-current absorption, pair effective opacity, and opt-in stiffness/projection thresholds need intermediate range beyond that of a materialized product. Published endpoints remain finite checked doubles; the ordinary source path preserves its existing rounding and zero-underflow policy. | `unit_test_m1_neutrino_source_update.c` checks scaled emission/absorption endpoints, `test_pair_scaled_opacity_fallback`, public policy selection, and helper boundary arms; `unit_test_m1_fd_jacobian.c` exercises endpoint projection recovery/rejection. |
| Normalized metric norms in `ghl_m1_utils.h` | Realizability and Eulerian/transport velocity require a norm even when the unnormalized quadratic form overflows or underflows. Normalized Cholesky evaluation and double accuracy gates avoid platform-specific acceptance. | `unit_test_m1_error_handling.c` checks private norm boundaries and realizability repair; `unit_test_m1_closure_fallback.c` and seeded invariants check nonzero fluxes and scaled states. |

Rejecting every nonfinite intermediate would change these established finite
endpoint contracts, not merely simplify their implementation. No blanket
range-extension promise is made: other overflowing contractions, ill-conditioned
metrics outside accuracy gates, and nonrepresentable final fluxes/endpoints still
return the documented errors. The existing `-ffp-contract=off` build route
checks these explicit-double paths without relying on implicit contraction.

## Build and run

From the repository root, for a build without HDF5:

```sh
./configure --noomp --disable-hdf5 --prefix="$PWD"
make test/unit_test_m1_closure_fallback test/unit_test_m1_diffusion_flux test/unit_test_m1_error_handling test/unit_test_m1_fd_jacobian test/unit_test_m1_neutrino_rusanov_flux test/unit_test_m1_neutrino_seeded_invariants test/unit_test_m1_neutrino_source_update test/unit_test_m1_rate_provider test/unit_test_m1_thcm1_blended_rusanov test/unit_test_rusanov_flux
```

For HDF5 support, omit `--disable-hdf5` and configure its installation normally.
Use a fresh build when changing compiler or configuration. These named
targets build only the M1 owners and their required library/test support.
Run local checks from the repository root with `LD_LIBRARY_PATH` pointing to
`build/lib`. For retained replay, pass `--fixture-dir PATH` to the seeded,
neutrino Rusanov, blended transport and generic Rusanov executables, and
`--fixture PATH` to the diffusion executable. Requested missing or incomplete
fixtures are errors. No Einstein Toolkit build is required.
The published package stores the versioned binary fixtures as individual
`.bin.gz` files under `radiation/`. Gzip is a storage wrapper: preparation restores
the original `.bin` names and `M1THCMB1` contents before the existing reader
runs. The reader also accepts the original text layout for local fixtures and
parser rejection checks.

With `M1_FIXTURE_DIR` unset, the ordinary runner downloads every named
`.bin.gz` member from `raw.githubusercontent.com` at the full TestData commit
in `.github/radiation-testdata-ref` through its existing shared download
machinery, then prepares the members in a private directory under
`${TMPDIR:-/tmp}` registered with its existing EXIT cleanup before
preparation. Supplying `M1_FIXTURE_DIR` bypasses the downloads and replays
the supplied raw or gzip directory while preserving caller inputs: a complete
raw directory is used directly, and compressed members are prepared in a
runner-owned temporary directory while the caller-owned source is left
untouched. Missing members and failed decompression
stop preparation before the stored-reference tests run. Direct individual
executable invocations without fixture arguments remain local-only regardless
of `M1_FIXTURE_DIR`.

To replay an existing sibling TestData checkout without downloads, run from
the GRHayL root:

```sh
M1_FIXTURE_DIR="$(pwd)/../TestData/radiation" \
  bash .github/run_tests.sh m1 --noomp
```

This uses only the local package. No publication, revision pin, or remote
fixture download is performed; the runner prepares compressed members once and
removes only its temporary output on exit.

For standalone replay, create an empty destination and retain ownership of its
cleanup:

Any raw or gzip local package works as the fixture source, including a
sibling `../TestData/radiation` checkout.

```sh
set -e
fixture_source="$(pwd)/../TestData/radiation"
fixture_dir=$(mktemp -d "${TMPDIR:-/tmp}/grhayl-m1.XXXXXX")
trap 'rm -rf -- "$fixture_dir"' EXIT
bash .github/prepare_m1_fixtures.sh "$fixture_source" "$fixture_dir"
export LD_LIBRARY_PATH="$PWD/build/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
./test/unit_test_m1_neutrino_seeded_invariants --fixture-dir "$fixture_dir"
./test/unit_test_m1_diffusion_flux --fixture "$fixture_dir/jthick_thcm1.bin"
```

Keep the destination until all requested replay executables finish. The same
helper accepts raw members and copies them into the supplied
empty destination. Gzip preparation preserves every record and binary64 value;
it changes neither the fixture grammar nor the numerical comparisons. The
stored package has a user-required ceiling of 10 MiB, including its sidecars;
expanded fixtures still require approximately 169.61 MiB outside TestData.

The ordinary runner uses `build/lib` and passes `--generated-fixture` to the
rate-provider test: HDF5 builds create and remove a
deterministic test-local EOS table and exercise the table-backed provider.
Without HDF5, the same invocation runs the available table-free checks. The
provider executable retains its optional external EOS-table argument. This
table is provider-test input, not an independent reference output.

An in-table cell is not necessarily representable by the rate provider. The
random sampler independently identifies a required negative Fermi--Dirac
exponential tail that rounds to zero (neutrino equilibrium moments, or
electron/positron moments when the pair channel is enabled). Those cells must
return the microphysics error without changing caller-owned rates or cache;
all other sampled cells retain the success and rate-bundle assertions.
Unexplained errors are not skipped or accepted. A successful table-backed
control and cached-tail rejection witnesses also run with the generated table.
Malformed-table witnesses modify and restore the queried midpoint's actual
interpolation corners, rather than assuming it lies in the first table cell.

## NRPyLeakage adapter scope

The M1 NRPyLeakage backend uses shared NRPyLeakage-owned inline helpers for
local nucleon blocking, reaction-shifted beta moments, and bremsstrahlung moments:

- `GRHayL/include/ghl_nrpyleakage_nucleon_blocking.h`
- `GRHayL/include/ghl_nrpyleakage_rate_helpers.h`

`GRHayL/include/make.code.defn` installs these headers. The adapter includes
this public boundary; helper implementations remain owned by NRPyLeakage.
M1 remains tau-free: optical-depth suppression, leakage luminosity/source
assembly, and leakage finite-output fallback policies do not enter the M1
rate bundle.

## Stored reference boundary

Explicit replay needs the payloads, not live THC_M1 execution. The fresh
candidate now includes producer records for every family. Publication in
`GRHayL/TestData` supplies the matched set at the full commit selected by
`.github/radiation-testdata-ref`; the
[TestData provenance](https://github.com/GRHayL/TestData/blob/main/radiation/PROVENANCE.md)
records the published package details. The ordinary runner downloads the
pinned `radiation/*.bin.gz` members by default, and the CI action forwards
its optional `fixture-dir` input to the runner as `M1_FIXTURE_DIR`; a nonempty
`fixture-dir` input instead replays a supplied local package without
downloads.
Current GRHayL outputs never become expected values. The staged receipt,
command, sidecar, consumed-input, and raw-result records are the evidence for
the new candidate; mutable offline campaign results are not a substitute for
those exact captured operands.

The exact variable-volume transport family is loaded by the existing
`unit_test_m1_thcm1_blended_rusanov` executable from
`transport_four_point_varying_d0.bin`, `transport_four_point_varying_d1.bin`,
`transport_four_point_varying_d2.bin`, and the
`transport_four_point_varying_controls.bin` control shard. It uses the existing 50-field
transport input layout: face volume, speed, opacity, spacing, theta, minimum
dissipation, four cell volumes, four five-component stencil states, and two
five-component physical face fluxes. The test adapter validates the volumes
and prepares `V_cell U` and `V_face F` before calling the required
`ghl_m1_compute_neutrino_four_point_volume_weighted_transport_flux`
transport API. For requested replay, missing variable-volume shards are a hard test failure; no
retained reference output is generated or substituted at runtime.

Reference comparisons supplement local tests. Matching a pointwise or prepared
face output does not establish full-grid evolution, continuum convergence, or
whole-code equivalence. Operations with different numerical or physical models
retain local tests and an explicit comparison limitation. The discrete-operator
equivalence workflow, live THC_M1 execution, host evolution, and provider/rate
workflows remain separate claims.

With explicit fixture arguments, these executables replay stored
baseline/perturbed outputs attributed to THC_M1:
`unit_test_m1_neutrino_seeded_invariants`,
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
evaluates current GRHayL at both inputs, preserves the retained baseline/response
envelope check, and checks both endpoints and the propagated paired response.
The propagated bound measures endpoint/response consistency, not independent
derivative accuracy. Pointwise,
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

The seeded-invariants owner also replays `stress_energy.bin`. Its complete
21-value state/ADM/velocity input is evaluated through the public neutrino
stress-energy wrapper at both paired endpoints, lowered with the full ADM
four-metric, and compared with the retained THC_M1 `assemble_rT` endpoints and
response. This remains a discrete stress-energy check; it does not claim host
evolution or continuum agreement.

There is no stored THC finite-step pair/implicit endpoint corpus. Independent
local collision-model oracles and solver checks cover those operations. The
explicit source-update, repair, and failure-path checks therefore remain local
evidence. See
[RUSANOV.md](m1_thcm1/RUSANOV.md) for the supplemental nonzero-current
fixture family, and [STRESS_ENERGY.md](m1_thcm1/STRESS_ENERGY.md) for the
covariant stress-energy corpus. Each consumer validates its requested fixture
grammar and operation-specific contracts before comparison.

## CI

The single consolidated [test runner](../../.github/run_tests.sh) takes
`[all|m1] [configure arguments...]`. `all` (the default) builds and runs every
unit test, including the M1 executables. `m1` builds with `make tests` and
runs only the Radiation M1 tests; the remaining arguments are passed to
`./configure -r`. With `M1_FIXTURE_DIR` unset, the runner refreshes the named
`radiation/*.bin.gz` members from `raw.githubusercontent.com` at the pinned
`.github/radiation-testdata-ref` revision through its existing shared download
machinery and runs the full stored replay by default. Supplying
`M1_FIXTURE_DIR` selects the supplied raw or gzip package and bypasses
downloads; direct individual executable invocations without fixture arguments
remain local-only. No Einstein Toolkit build is required.

Each compiler workflow in
[`.github/workflows/`](../../.github/workflows/) (Ubuntu gcc, clang and Intel;
macOS gcc and clang) has a `radiation-m1` job with an OS and HDF5
enabled/disabled matrix. The job calls the
[`run_m1` action](../../.github/actions/run_m1/action.yml), which selects the
compiler, the HDF5 mode (`--disable-hdf5` for the disabled leg), the Homebrew
compilers on macOS, and gcc coverage flags on Linux gcc, passes `--noomp`, and
ends with `.github/run_tests.sh m1`. Only the
[Ubuntu GCC workflow](../../.github/workflows/github-actions-Ubuntu-gcc.yml)
job adds a coverage step: the
[M1 coverage action](../../.github/actions/m1-code-coverage/action.yml)
generates a gcovr Cobertura report filtered to the Radiation sources plus
`ghl_m1.h` and `abort_if_error.c`; the Rusanov flux source is already under
Radiation. It uploads only `m1-coverage.xml` under the `m1` flag, disables
upload discovery, and fails on upload errors. The Radiation Codecov component
selects that flag so unrelated zero-hit reports cannot dilute its coverage.
The other `radiation-m1` jobs are plain test runs without a coverage upload.
Legacy jobs retain the shared coverage action's existing collection and upload
discovery, without M1-specific exclusions or gates.

The action's Radiation gate uses the `GRHayL/Radiation/` source filter,
recursively including `Neutrinos/` and executable private headers, and
requires 100% executable lines; branch coverage is reported without a gate.
Branch coverage resolves partial line coverage; this is not MC/DC or path
coverage. Functions are reported as supporting evidence without another gate.
Public inline helpers and Core error mappings are supplemental evidence and do
not inflate the Radiation denominator. An excluded arm is not an executed arm:
a passing gate that relies on an exclusion requires an independently checked
reachability or compiler-accounting proof for that exemption.

The fixture-integration contract for the `radiation-m1` CI jobs is
owned by the CI files: the `run_m1` action replays the
pinned TestData publication by default through the existing runner route:
it forwards its optional `fixture-dir` input to the runner as
`M1_FIXTURE_DIR`, and the runner performs the pinned
`.github/radiation-testdata-ref` acquisition. The revision file
records a full TestData commit, not an occurrence checksum. A nonempty
`fixture-dir` input instead replays the supplied raw or gzip package and
skips downloads. External acquisition fails if the pinned commit is
unavailable in `GRHayL/TestData`; the action does not fall
back to a moving branch or skip stored replay. Direct individual executable
invocations without fixture arguments remain local-only. Workflow
configuration and local runner verification do not establish that a remote
CI run has passed.

For local action replay, supply its optional nonempty `fixture-dir` input
with an existing raw or gzip fixture directory. The default empty input
selects the pinned publication route, whose revision file must then contain
a full published TestData commit SHA; a missing or malformed pin fails with
a publication-prerequisite diagnostic and never silently omits stored replay.

The provider tests distinguish test-local reference implementations from real
production calls. Enabled builds exercise the real table-backed provider with
`--generated-fixture`; disabled builds exercise its supported error boundary
and the private disabled-backend stubs, including unchanged outputs and
validation precedence. Table/cache/recovery operations are compiled only in
the supported enabled implementation. The transport test resets the actual
stencil row extent; its focused UBSan evidence uses counters separate from the
coverage measurement.

The source-update fault-injection cases in
`Unit_Tests/m1_source_update_fault_injection.h` and the repair post-norm cases
in `Unit_Tests/unit_test_m1_error_handling.c` use Linux ELF interposition and
`dlsym(RTLD_NEXT)`. Untargeted calls delegate to the real library. Their callers
guard them with `__linux__`; macOS builds execute the ordinary checks without
these interposition cases. The FD-Jacobian owner also exercises pair quadratic
range rejection, scaled endpoint-number recovery, and the generic Newton
driver's own admissibility and metric checks after successful callbacks.
M1 exposes no public hook pointers for these cases: mutable/global hooks would
conflict with the caller-owned state in the public API. ELF interposition is
the accepted Linux-only injection mechanism; the ordinary source-update
coverage remains portable to macOS.

## Separate workflows and evidence limits

The following remain outside this unit-test claim boundary:

- live THC_M1 execution and mutable offline campaign results;
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

## Replay comparison policies

The normalized trusted/perturbed bar delegates to
`ghl_pert_test_fail_with_tolerance`. Its adapter preserves the established
absolute-bound equality, zero-reference and floating-point normalization
boundary behavior where the shared helper differs.
Zero-reference values are normalized to the recorded input scale before the
canonical helper delegation, so the helper's absolute zero-reference semantics
reproduce the established normalized response rule.

The endpoint/response campaign comparator remains separate: it is symmetric in
current/reference endpoint magnitudes, uses the retained `1e-300` magnitude
floor, and propagates endpoint uncertainty into the paired-response bound.
The shared trusted-baseline helper is asymmetric and cannot express that rule.
These are existing fixture policies; no tolerance is changed. The policy
strings (such as `source_a1_a2_rate_normalized_v1` and
`strict_relative_2e-12_propagated_response_v1`) and their test-side tolerance
constants are retained campaign metadata intentionally mirrored from the
recorded producer policy. They document the stored fixture's policy identity
and are not production tolerances.
