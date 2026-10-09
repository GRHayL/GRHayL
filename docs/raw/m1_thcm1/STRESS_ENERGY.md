# THC_M1 stress-energy fixtures

The prepared TestData package contains `stress_energy.bin` with 1,024 complete
baseline/perturbed pairs and ten covariant tensor components per pair. Its
historical notes attribute the values to the external CL-04 `stress_energy`
route and compiled THC_M1 `assemble_rT` outputs. The package manifest marks
them historical and not admitted for current verification; the receipt, source
snapshot, consumed inputs, and raw outputs needed to verify that attribution
were unavailable in the reviewed workspace. The current test recomputes GRHayL
at both endpoints and compares the baseline, perturbed endpoint, and
current-vs-current response with the retained values.

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

This payload is not present in the current GRHayL checkout or pinned to a
published TestData revision. The `.github/run_tests.sh` runner does not fetch
M1 data: it runs local-only M1 checks by default, and an explicitly supplied
`M1_FIXTURE_DIR` requests replay from local files. The retained sibling package
manifest declares its digest purpose as fixture-payload integrity only and
marks the package historical and unadmitted; those fields do not authenticate
the producer or the payload's campaign binding.
