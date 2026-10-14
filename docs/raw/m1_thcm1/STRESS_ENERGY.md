# THC_M1 stress-energy fixtures

The published TestData package contains `stress_energy.bin` with 1,024 complete
baseline/perturbed pairs and ten covariant tensor components per pair. The
values come from the fresh CL-04 `stress_energy` campaign, which evaluates the
compiled THC_M1 `assemble_rT` route over frozen baseline and perturbed packets
and reproduces every captured THC_M1 value exactly. The campaign's receipt,
source snapshot, consumed inputs, and raw outputs are in the producer evidence
outside this checkout; its original comparison metadata and pass flags are
retained rather than recomputed and are not new comparison evidence. The current
test recomputes GRHayL at both endpoints and compares the baseline, perturbed
endpoint, and current-vs-current response with the retained values.

Each consumed input has 21 values in this order:

```text
N, E, F[0], F[1], F[2], lapse, shiftU[0..2], gammaDD[0..8], vU[0..2]
```

`gammaDD` is row-major. The output order is
`T_dd[0][0]`, `T_dd[0][1]`, `T_dd[0][2]`, `T_dd[0][3]`, `T_dd[1][1]`,
`T_dd[1][2]`, `T_dd[1][3]`, `T_dd[2][2]`, `T_dd[2][3]`, and `T_dd[3][3]`.
GRHayL constructs the contravariant radiation tensor from E/F and the closure,
then the test lowers it with the full ADM four-metric before comparison.

Every pair changes only the radiation energy input, records that sensitivity
as `radiation.E`, and recomputes its positive normalization from the paired
input energies. The fixture uses the
`strict_relative_2e-12_propagated_response_v1` policy shared with the
prepared-transport fixtures: each current baseline and perturbed output must
agree with the retained THC value to relative 2e-12, and the response bound
propagates those endpoint bounds. It is a discrete stress-energy operation
check, not a full evolution, host-integration, or continuum-equivalence claim.

This payload is not present in the GRHayL checkout. The `.github/run_tests.sh m1`
runner downloads `stress_energy.bin.gz` and the other M1 members from the
TestData revision pinned in `.github/radiation-testdata-ref`, decompresses them
into a private temporary directory, and replays them. Setting `M1_FIXTURE_DIR`
to a local raw or gzip package skips the download, and running the executable
without fixture arguments runs local-only checks. The producer evidence's package
manifest declares its digest purpose as fixture-payload integrity only; the
digests do not authenticate the producer.
