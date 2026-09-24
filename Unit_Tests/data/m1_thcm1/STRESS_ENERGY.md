# THC_M1 stress-energy fixtures

`stress_energy.m1` is a portable, frozen replay of the external CL-04
`stress_energy` campaign. It contains 1,024 complete baseline/perturbed pairs
and ten covariant tensor components per pair. The retained values are the
compiled THC_M1 `assemble_rT` outputs; the current public test recomputes
GRHayL at both endpoints and compares the baseline, perturbed endpoint, and
current-vs-current response with the retained THC_M1 values.

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

The campaign producer/exporter lives in the external Verification workspace;
the public repository stores only this portable payload and compact manifest
metadata. Normal unit-test execution does not require THC_M1, `THCM1_ROOT`,
Verification, or a reference-data download. The package manifest intentionally
retains its existing historical/unadmitted package-level status because the
legacy payloads remain in that same package.
