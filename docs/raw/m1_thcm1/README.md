# Frozen Radiation THC_M1 reference data

Existing Radiation tests replay explicitly supplied fixtures without THC_M1,
`THCM1_ROOT`, or Verification execution. Offline exporters are not invoked by
normal tests. Current GRHayL outputs never become trusted expected values.

A fresh replacement candidate was generated for every family from new THC_M1
runs. Its manifest includes Jthick and its producer evidence
contains receipts, captured commands, consumed inputs, raw outputs, source
snapshots, and exporters. The published fixture families and binary layout
are described by the maintained
[TestData provenance](https://github.com/GRHayL/TestData/blob/main/radiation/PROVENANCE.md).
The accepted payloads are published as individual `radiation/*.bin.gz` members
in `GRHayL/TestData`; the ordinary runner downloads and prepares every member
at the pinned revision by default, the CI action retries that acquisition and
forwards its optional `fixture-dir` input to the runner as `M1_FIXTURE_DIR`.
A supplied local raw or gzip package replays without downloads or publication.
The [production record](../Radiation_fresh_fixture_production.md)
summarizes the scoped runs and their limits. The
producer-evidence archive stays outside this checkout.
The historical attribution discussion below describes the removed historical
package; it is not the provenance claim for the accepted one.

## Pointwise closure, moments, stress, geometry, and speeds

`pointwise_closure_moments.bin` is consumed by
`unit_test_m1_neutrino_seeded_invariants.c` alongside its existing local checks.
It contains 56 retained pairs. Historical family notes attribute them to an
external THC_M1 campaign capture and describe receipt, command-output, and
consumed-input binding required for admission. That earlier external pointwise
campaign retains its packet stream, THC_M1-labeled raw results, commands,
receipt, and source snapshot in its own workspace. All
112 retained input and output vectors match that campaign exactly. No exporter
receipt binds the campaign to this fixture package, so the historical producer
attribution remains unadmitted. The package manifest marks the payload
unadmitted; its retained values are not current-GRHayL-generated expectations.

The 68 input scalars are E, covariant F[3], lapse, shift[3], gammaDD[9],
coordinate velocity[3], K[9], lapse derivatives[3], shift derivatives[9], and
metric derivatives[27], with matrices in row-major order. The retained notes
record M1 controls `(1e-8,1e-12,1,1e-6,1e-10,100,1e-10)`; replay uses those same
controls and independently computes GRHayL's closure for downstream operations.
The 38 output scalars are P[9], chi, xi, J, covariant H[3], contravariant T4[16],
densitized geometry sources[4], and maximum light-cone speeds[3].

The archived `radiation/*.bin.gz` members expand to the versioned binary layout
documented in the maintained
[TestData provenance](https://github.com/GRHayL/TestData/blob/main/radiation/PROVENANCE.md):
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

The historical notes state that the producer stream also carries N. N is
omitted from this operation fixture because these APIs do not consume it; the
producer stream needed to verify its original binding is unavailable.
Twenty-nine pairs exercise actual pointwise input perturbations. Twenty-five change
only fields outside this operation's inputs and are labeled
`input_invariant_control`; they establish fixed-input agreement only.

Two named pairs have a documented policy difference:
`rngpkt-v2-rd-radiation-anchor-energy-floor-a01` and
`rngpkt-v2-rd-metric-anchor-offdiagonal-spd-a01`. Both have exactly zero Eulerian
flux and moving fluid. Their THC pressure fails the pressure-trace requirement
in `GRHayL/Radiation/ghl_m1_utils.h`. Current GRHayL deliberately publishes an
Eulerian Minerbo admissibility fallback. The test names these exceptions,
checks the reference trace discrepancy, and asserts current fallback status,
finite outputs, pressure trace, and isotropic closure factor. They count as local
admissibility checks, never THC agreement; they retain their explicit baseline
and perturbed local checks. A new or changed classification fails.

## Stress-energy tensor

`stress_energy.bin` is consumed by the existing
`unit_test_m1_neutrino_seeded_invariants` owner. It contains 1,024 complete
baseline/perturbed pairs and the ten symmetric covariant components. The
retained sibling package labels the operation `stress_energy`; its accompanying
notes attribute the records to the external CL-04 route and compiled THC_M1
`assemble_rT` outputs. The receipt, source snapshot, consumed inputs, and raw
producer output needed to verify that binding were unavailable in the reviewed
workspace. The package remains historical and unadmitted. The 21 public input
values are N, E, covariant F[3], lapse,
shift[3], row-major gammaDD[9], and coordinate fluid velocity[3]. Every pair
perturbs only E. See [STRESS_ENERGY.md](STRESS_ENERGY.md) for the component
mapping and evidence boundary.

The replay recomputes GRHayL's stress-energy tensor at both paired inputs and
lowers each result with the full ADM four-metric. It compares both endpoints
and the current-vs-current response with retained endpoints and responses
attributed to THC_M1 `assemble_rT`; all 1,024 paired records match under the
stored policy. The missing producer binding prevents treating that match as
an admitted THC_M1 comparison. THC_M1, `THCM1_ROOT`, and offline producer
campaigns remain offline provenance inputs rather than public test
dependencies.

## Other operations

See [SOURCE.md](SOURCE.md) for the frozen instantaneous-source inputs, 162
agreement pairs, and six named local-policy cases. See
[TRANSPORT.md](TRANSPORT.md) for prepared-transport provenance, retained
input conventions, and the strict paired endpoint/response comparison policy.
See [RUSANOV.md](RUSANOV.md) for the two complete 1024-face paired corpora
and their existing generic/neutrino Rusanov test owners.
See [STRESS_ENERGY.md](STRESS_ENERGY.md) for the complete covariant
stress-energy corpus and its public tensor mapping.

Neither the removed historical package nor the accepted payloads live in this
checkout. The candidate's producer-evidence archive stays outside the source
tree with its family manifest, including Jthick, and its
`fixture/README.md` maps each family to its fresh producer records and
exporter. The maintained
[TestData provenance](https://github.com/GRHayL/TestData/blob/main/radiation/PROVENANCE.md)
is the published package authority; the old package remains historical and is
not the source of the new producer attribution.

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

The retained family notes describe pointwise results emitted in consumed-input
order, without row IDs. Any admission must bind that order to the actual
recorded input stream. Copying a workspace package or replaying current GRHayL
does not establish that historical binding.
