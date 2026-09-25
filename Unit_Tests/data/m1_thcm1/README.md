# Frozen Radiation THC_M1 reference data

Existing Radiation tests consume these fixtures without THC_M1, `THCM1_ROOT`,
Verification, or reference-data downloads. Offline exporters are not invoked by
normal tests. Current GRHayL outputs never become trusted expected values.

## Pointwise closure, moments, stress, geometry, and speeds

`pointwise_closure_moments.m1` is consumed by
`unit_test_m1_neutrino_seeded_invariants.c` alongside its existing local checks.
It preserves 56 pairs from a retained external THC_M1 campaign capture. The
external exporter binds the campaign receipt, recorded command output, and
consumed-input stream before any future promotion. Those source artifacts are
not part of the public GRHayL repository, and the retained provenance does not
assert that current GRHayL is identical to the historical producer snapshot.

The 68 input scalars are E, covariant F[3], lapse, shift[3], gammaDD[9],
coordinate velocity[3], K[9], lapse derivatives[3], shift derivatives[9], and
metric derivatives[27], with matrices in row-major order. The frozen producer
initializes M1 with `(1e-8,1e-12,1,1e-6,1e-10,100,1e-10)`; replay uses those same
controls and independently computes GRHayL's closure for downstream operations.
The 38 output scalars are P[9], chi, xi, J, covariant H[3], contravariant T4[16],
densitized geometry sources[4], and maximum light-cone speeds[3].

The versioned token format starts with `M1_THCM1_FIXTURE 1`, operation, policy,
and record count. Each named record stores IDs and perturbation metadata, both
input vectors, both THC output vectors, normalization, and availability/status
fields; `record_end` and final `end` terminate the sections. The test-local
`m1_thcm1_fixture_utils.h` defines the exact grammar and validates counts,
finiteness, IDs, and successful reference availability. Float text uses round-trip
precision. Consumers independently recompute normalization from the inputs.

The `pointwise_a1_a2_v1` policy is the retained verification policy: absolute
`2e-12`, relative `2e-10`, floor `1e-13` after input-defined normalization.
Pressure, moments, stress, and geometry use the larger input E; closure scalars
and speeds use unit normalization. The ordinary pointwise replay evaluates
current GRHayL at both baseline and perturbed inputs, and compares its endpoint
response directly with the retained THC response. These are reference-test
bounds, not production solver settings.

The producer stream also carries N; its bytes are retained as part of the
recorded input binding, but N is omitted from the operation fixture because
these APIs do not consume it.
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

`stress_energy.m1` is consumed by the existing
`unit_test_m1_neutrino_seeded_invariants` owner. It contains 1,024 complete
baseline/perturbed pairs from the external CL-04 `stress_energy` route and
the ten symmetric covariant components emitted by compiled THC_M1
`assemble_rT`. The 21 public input values are N, E, covariant F[3], lapse,
shift[3], row-major gammaDD[9], and coordinate fluid velocity[3]. Every pair
perturbs only E. See [STRESS_ENERGY.md](STRESS_ENERGY.md) for the component
mapping and evidence boundary.

The replay recomputes GRHayL's stress-energy tensor at both paired inputs and
lowers each result with the full ADM four-metric. It compares both endpoints
and the current-vs-current response with the retained THC_M1 `assemble_rT`
endpoints and response; all 1,024 paired stress-energy responses agree. THC_M1,
`THCM1_ROOT`, and `Verification/` remain
offline provenance inputs rather than public test dependencies.

## Other operations

See [SOURCE.md](SOURCE.md) for the frozen instantaneous-source inputs, 162
agreement pairs, and six named local-policy cases. See
[TRANSPORT.md](TRANSPORT.md) for prepared-transport provenance, admitted
input conventions, and the strict paired endpoint/response comparison policy.
See [RUSANOV.md](RUSANOV.md) for the two complete 1024-face paired corpora
and their existing generic/neutrino Rusanov test owners.
See [STRESS_ENERGY.md](STRESS_ENERGY.md) for the complete covariant
stress-energy corpus and its public tensor mapping.

The tracked [payload archive](payloads.tar.gz) contains the 13 plaintext
fixtures listed in [package_manifest.json](package_manifest.json). A clean
checkout's runner expands it once, then validates package membership, records,
and payload digests with [audit_package.py](audit_package.py) before replay.
The expanded `.m1` and `.dat` files are ignored by git. Existing plaintext
payloads are audited without replacement.

The pointwise producer emits results in consumed-input order, without row IDs.
The external promotion step binds that order to the recorded input stream.
Public package digests detect payload mutation; they are not source tracking or
cryptographic authentication of a historical producer. The current package
manifest remains explicit about its historical/unadmitted status until an
external admission receipt is supplied.
