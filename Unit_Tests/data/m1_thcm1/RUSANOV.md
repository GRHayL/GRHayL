# Paired Rusanov fixtures

These files contain the complete retained 1024-face baseline/perturbed corpus
from the `rusanov_efn_flux_low_matched` producer route:

- `rusanov_neutrino.dat`: `neutrino_rusanov_flux`, 54 consumed inputs and five
  `{N,E,Fx,Fy,Fz}` outputs.
- `rusanov_generic.dat`: `rusanov_flux`, 17 consumed inputs and four `{E,Fx,Fy,Fz}`
  outputs.
- `rusanov_neutrino_current.dat`: `neutrino_rusanov_flux_current_v2`, the same
  54-input/five-output layout with nonzero number-current operands in slots
  41--52.  Slots 41--43 and 44--46 are `N*velocity` for the left and right
  states; slots 47--49 and 50--52 retain the corresponding transport
  velocities.

The two original files use `strict_relative_2e-12_propagated_response_v1`, retain every face
pair, and require face IDs in order `face0000` through `face1023`.  The
retained direction coverage is 342 faces in direction 0, followed by 341 in
directions 1 and 2.

The current shard uses
`strict_relative_2e-12_propagated_response_current_v2` and retains all 1024
pairs under the separate operation identity
`neutrino_rusanov_flux_current_v2`.  Its `radiation_N` perturbations change
the stored THC number output, so removing, reversing, or rescaling the
current operands is observable during fixture evaluation.

All 1024 supplemental pairs change consumed inputs and outputs; 341 are
number-density perturbations. The original neutrino corpus has 1024 changed
pairs, while the generic projection has 683 changed pairs and 341 controls.

The retained corpora use the current public GRHayL Rusanov APIs together with
THC_M1 face-local operands and outputs. The fixture stores caller-prepared
metric, closure, number-current, velocity, speed, and physical-flux operands
as one serialized common input. It therefore exercises the Rusanov operation
against the stored baseline and perturbation response; it is not independent
validation of closure or physical-flux preparation. The generic fixture
likewise stores its caller-supplied physical E/F operands.

These payloads are currently retained historical data, not a current admitted
cross-code result. External campaign generation, admission, receipts, raw
captures, and source/build context remain outside the public GRHayL checkout.
When a payload is admitted, only its portable fixture data and compact
receipt/artifact reference are promoted into this directory. Normal Unit_Tests
never invoke the external campaign or require its checkout.

`Unit_Tests/m1_thcm1_rusanov_fixture.h` exposes the test-local consumers:

```c
m1_thcm1_rusanov_check_neutrino_fixture(...)
m1_thcm1_rusanov_check_current_fixture(...)
m1_thcm1_rusanov_check_generic_fixture(...)
```

The two existing Rusanov unit tests invoke these consumers before their
pre-existing local checks and accept `--fixture-dir PATH`.  THC_M1 is not a
runtime dependency of either unit test.
