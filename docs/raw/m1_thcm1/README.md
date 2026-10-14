# Frozen Radiation THC_M1 reference data

Existing Radiation tests replay explicitly supplied fixtures without THC_M1,
`THCM1_ROOT`, or Verification execution. Offline exporters are not invoked by
normal tests. Current GRHayL outputs never become trusted expected values.

Every family was produced from fresh THC_M1 runs, including the standalone
Jthick scalar. The producer evidence contains receipts, captured commands,
consumed inputs, raw outputs, source snapshots, and exporters. The published
fixture families, producer, regeneration route, and binary layout are described
by `radiation/PROVENANCE.md` in `GRHayL/TestData` at the revision pinned in
`.github/radiation-testdata-ref`.
The payloads are published as individual `radiation/*.bin.gz` members
in `GRHayL/TestData`; the ordinary runner downloads and prepares every member
at the pinned revision by default, the CI action retries that acquisition and
forwards its optional `fixture-dir` input to the runner as `M1_FIXTURE_DIR`.
A supplied local raw or gzip package replays without downloads or publication.
The [production record](../Radiation_fresh_fixture_production.md)
summarizes the scoped runs and their limits. The
producer-evidence archive stays outside this checkout. An earlier historical
package that lacked a verifiable producer-to-payload binding has been removed;
see [Superseded historical package](#superseded-historical-package).

## Pointwise closure, moments, stress, geometry, and speeds

`pointwise_closure_moments.bin` is consumed by
`unit_test_m1_neutrino_seeded_invariants.c` alongside its existing local checks.
It contains 56 retained pairs. They come from the fresh pointwise THC_M1
campaign, whose packet stream, consumed inputs, commands, receipt, raw results,
and source snapshot are in the producer evidence outside this checkout. The
retained values are THC_M1 outputs, not current-GRHayL-generated expectations.

The 68 input scalars are E, covariant F[3], lapse, shift[3], gammaDD[9],
coordinate velocity[3], K[9], lapse derivatives[3], shift derivatives[9], and
metric derivatives[27], with matrices in row-major order. Replay initializes
GRHayL with the M1 controls `(1e-8,1e-12,1,1e-6,1e-10,100,1e-10)` and
independently computes GRHayL's closure for downstream operations.
The 38 output scalars are P[9], chi, xi, J, covariant H[3], contravariant T4[16],
densitized geometry sources[4], and maximum light-cone speeds[3].

The archived `radiation/*.bin.gz` members expand to the versioned binary layout
documented in the pinned `radiation/PROVENANCE.md` of `GRHayL/TestData`:
eight-byte magic `M1THCMB1`, little-endian
integer and binary64 fields, and length-prefixed UTF-8 strings and vectors.
The original text export starts with `M1_THCM1_FIXTURE 1`, operation, policy,
and record count. Each named record stores IDs and perturbation metadata, both
input vectors, both THC output vectors, normalization, and availability/status
fields; `record_end` and final `end` terminate the sections. The test-local
`m1_thcm1_fixture_utils.h` validates both layouts, counts, finiteness, IDs,
and successful reference availability. The binary files preserve the text
export's binary64 values. Consumers independently recompute normalization
from the inputs.

The `pointwise_a1_a2_v1` policy is the retained comparison policy: absolute
`2e-12`, relative `2e-10`, floor `1e-13` after input-defined normalization.
Pressure, moments, stress, and geometry use the larger input E; closure scalars
and speeds use unit normalization. The ordinary pointwise replay evaluates
current GRHayL at both baseline and perturbed inputs, and compares its endpoint
response directly with the retained THC response. These are reference-test
bounds, not production solver settings.

The pointwise producer stream carries one more scalar, N, than the 68 stored
inputs. N is omitted from this operation fixture because these APIs do not
consume it.
Thirty-one pairs exercise actual pointwise input perturbations: 28 ordinary
strict comparison pairs and the three named fallback-policy pairs described
below. The other 25 change only fields outside this operation's inputs and are
labeled `input_invariant_control`; they establish fixed-input agreement only.

Three named pairs have a documented policy difference, carried as independent
local-policy classifications rather than cross-code replays. Two,
`rngpkt-v2-rd-radiation-anchor-energy-floor-a01` and
`rngpkt-v2-rd-metric-anchor-offdiagonal-spd-a01`, have exactly zero Eulerian
flux and moving fluid. Their THC pressure fails the pressure-trace requirement
in `GRHayL/Radiation/ghl_m1_utils.h`, and current GRHayL deliberately publishes
an exact-zero-flux Eulerian Minerbo admissibility fallback. The third,
`rngpkt-v2-rd-radiation-anchor-realizability-cone-a01`, retains a THC pointwise
record whose baseline and perturbed endpoints share the identity spatial metric
and an x-aligned flux but whose reference pressure is nonrealizable: the
Eulerian covariance C^{ij} = P^{ij} - F^iF^j/E of that reference is negative
along the radiation flux direction. This is a fallback-policy distinction,
separate from the exact-zero-flux cause above and from any refreshed trusted
expected data: current GRHayL publishes its finite Eulerian covariance
admissibility fallback instead of reproducing the nonrealizable reference
pressure. The test names
each exception: the zero-flux records keep their established checks of
fallback status, finite outputs, reference pressure trace, and isotropic
closure factor. The covariance record is classified by a negative-Cxx/E
witness computed directly from the retained reference, independently of the
production covariance validator; for it the test separately validates current
fallback status and metadata (zero root iterations and residual), a finite
symmetric current fallback pressure with a valid trace, diagonal covariance
positivity, and the analytic Eulerian Minerbo fallback pressure and closure
factor. These three records are therefore not numeric cross-code comparisons:
their stored THC outputs serve only as the classification evidence named above,
and a small change to them need not fail the test. All other pointwise records
retain strict cross-code paired replay. A new or changed classification fails.

## Stress-energy tensor

`stress_energy.bin` is consumed by the existing
`unit_test_m1_neutrino_seeded_invariants` owner. It contains 1,024 complete
baseline/perturbed pairs and the ten symmetric covariant components. The
records come from the fresh CL-04 stress-energy campaign, which evaluates the
compiled THC_M1 `assemble_rT` route over frozen baseline and perturbed packets;
its receipt, source snapshot, consumed inputs, and raw producer output are in
the producer evidence outside this checkout. The 21 public input
values are N, E, covariant F[3], lapse,
shift[3], row-major gammaDD[9], and coordinate fluid velocity[3]. Every pair
perturbs only E. See [STRESS_ENERGY.md](STRESS_ENERGY.md) for the component
mapping and evidence boundary.

The replay recomputes GRHayL's stress-energy tensor at both paired inputs and
lowers each result with the full ADM four-metric. It compares both endpoints
and the current-vs-current response with retained endpoints and responses
produced by THC_M1 `assemble_rT`; all 1,024 paired records match under the
stored policy. The match is a discrete-operation comparison against stored
THC_M1 values, not a host-evolution claim. THC_M1, `THCM1_ROOT`, and offline
producer campaigns remain offline provenance inputs rather than public test
dependencies.

## Other operations

See [SOURCE.md](SOURCE.md) for the frozen instantaneous-source inputs, 159
agreement pairs, six named exact-zero-flux local-policy cases, and the named
covariance-anchor derived-species policy cases. See
[TRANSPORT.md](TRANSPORT.md) for prepared-transport provenance, retained
input conventions, and the strict paired endpoint/response comparison policy.
See [RUSANOV.md](RUSANOV.md) for the three complete 1024-face paired corpora
and their existing generic/neutrino Rusanov test owners.
See [STRESS_ENERGY.md](STRESS_ENERGY.md) for the complete covariant
stress-energy corpus and its public tensor mapping.

`jthick_thcm1.bin` holds six baseline/perturbed pairs of the thick-limit scalar
`Jthick` (20 input scalars, one output, policy `pointwise_a1_a2_v1`, each pair
perturbing `radiation.E`). It is consumed by
`unit_test_m1_diffusion_flux --fixture`. The Jthick helper campaign receipt
records 55 agreements and one difference between the campaign's GRHayL build and
THC_M1; that difference is a recorded limitation of the campaign comparison,
and the fixture stores only THC_M1 outputs. The optional Fick flux remains
local-only.

The accepted payloads do not live in this checkout. The producer-evidence
archive stays outside the source tree with its family manifest, including
Jthick, and its `fixture/README.md` maps each family to its fresh producer
records and exporter. The pinned `radiation/PROVENANCE.md` of `GRHayL/TestData`
is the published package reference: it lists the fixture families, producer,
regeneration route, and binary layout.

Replay requires every named member, including the standalone Jthick fixture.
With `M1_FIXTURE_DIR` unset, the ordinary `.github/run_tests.sh m1` runner
refreshes every named `.bin.gz` member from the pinned
`.github/radiation-testdata-ref` revision of `GRHayL/TestData` through the
existing download machinery and prepares the members in a private
`${TMPDIR:-/tmp}` directory removed by its EXIT cleanup; direct individual
executable invocations without fixture arguments remain local-only. Setting
`M1_FIXTURE_DIR` to a supplied raw or gzip package bypasses downloads and
preserves caller inputs, so a sibling local TestData checkout replays locally
without publication or a revision pin. The CI action retries the pinned
acquisition by default, forwards its optional nonempty `fixture-dir` input to
the runner as `M1_FIXTURE_DIR`, and selects the supplied package only when
that input is nonempty.

## Superseded historical package

An earlier package, whose records were attributed to external campaigns by
historical notes, has been removed. Its receipts, source snapshots, consumed
inputs, and raw producer output were unavailable, so its producer-to-payload
binding could never be verified, and replaying current GRHayL against it could
not establish that binding. None of its attribution claims is carried forward:
the provenance of the current payloads rests only on the fresh producer
evidence described above.
