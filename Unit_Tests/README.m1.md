# Radiation M1 unit tests

The two-layer evidence boundary for this directory is summarized in
[`README.md`](README.md). This page records the current scoped M1 inventory,
stored fixture families, and their evidence limits.

The existing M1 test executables exercise shared Radiation kernels and grey
three-species neutrino operations. Their randomized, analytic, invariant, and
transactional checks remain part of normal execution. Stored THC_M1 comparisons
are incorporated into the relevant existing executables; fixture provenance and
operation-specific coverage are documented with the data in `data/m1_thcm1/`.

Local regressions in `unit_test_m1_neutrino_seeded_invariants` cover nonzero
fluxes whose squared norms underflow or whose `E/F^2` intermediates overflow,
while retaining the distinct exact-zero closure policy. The
`unit_test_m1_neutrino_source_update` regressions check mean-energy bounds at
the final paired endpoint, including valid endpoints reached through an
out-of-bounds intermediate state and transactional rejection of invalid final
endpoints.

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

The runner uses `build/lib` for the current checkout's shared library. It passes
`--generated-fixture` to the rate-provider test: HDF5 builds create and remove a
deterministic test-local EOS table and exercise the table-backed provider.
Without HDF5, the same invocation runs the available table-free checks. The
provider executable retains its optional external EOS-table argument. This
table is provider-test input, not an independent reference output.

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
package; they do not re-run or replace Verification admission. This checkout
does not include a current Verification admission manifest or mutable campaign
results, so no current target admission, target alias, or fresh campaign status
is asserted by an M1 unit-test pass. Fixture refresh or promotion remains
blocked until Verification supplies its accepted receipt and status.

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
trusted output is generated or substituted at runtime.

Reference comparisons supplement local tests. Matching a pointwise or prepared
face output does not establish full-grid evolution, continuum convergence, or
whole-code equivalence. Operations with different numerical or physical models
retain local tests and an explicit comparison limitation. The discrete-operator
equivalence workflow, live THC_M1/Verification execution, host evolution, and
provider/rate workflows remain separate claims.

Four of the nine executables replay stored trusted baseline/perturbed THC
outputs: `unit_test_m1_neutrino_seeded_invariants`,
`unit_test_m1_neutrino_rusanov_flux`, `unit_test_m1_thcm1_blended_rusanov`, and
`unit_test_rusanov_flux`. The established corpus has these distinctions:

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

Pointwise and instantaneous-source families include respectively two and six
named exact-zero-flux moving-fluid policy differences. Those records exercise
local policy checks; they are not THC agreement claims. Identical-input records
are controls, and changed inputs need not change every output. Ordinary stored
replay evaluates current GRHayL once at each baseline input; the retained THC
perturbation supplies the response envelope and is not recomputed by GRHayL.
The prepared face operands also limit
Rusanov/transport agreement to the discrete operation, rather than independently
validating physical flux and closure preparation.

For the pointwise owner, ordinary records use one current GRHayL evaluation at
the baseline input. The paired THC_M1 perturbed output is the recorded response
envelope; the current GRHayL perturbed endpoint is not recomputed. The two named
zero-flux policy records remain explicit local two-state admissibility checks.

The seeded-invariants owner also replays `stress_energy.m1`. Its complete
21-value state/ADM/velocity input is evaluated through the public neutrino
stress-energy wrapper, lowered with the full ADM four-metric, and compared with
the retained THC_M1 `assemble_rT` baseline. The paired THC_M1 perturbation of
E supplies the strict local response envelope. This remains a discrete
stress-energy check; it does not claim host evolution or continuum agreement.

There is no stored THC finite-step pair/implicit endpoint corpus. Independent
local collision-model oracles and solver checks cover those operations. The
explicit source-update, repair, and failure-path checks therefore remain local
evidence. See
[RUSANOV.md](data/m1_thcm1/RUSANOV.md) for the supplemental nonzero-current
fixture family, and [STRESS_ENERGY.md](data/m1_thcm1/STRESS_ENERGY.md) for the
covariant stress-energy corpus. The runner validates the retained package with
`audit_package.py` before executing its consumers.

## CI and coverage

`.github/run_tests.sh` invokes the M1 runner. The existing compiler/OS workflows
also select dedicated Radiation jobs through `.github/actions/run_m1`; these
build the same named targets and execute both HDF5 configurations. Linux GCC jobs
retain a Radiation-filtered coverage JSON artifact and submit the executed M1
coverage through the existing Codecov action to the PR patch gate. Workflow selection is not a
claim that a remote CI run has already passed.

Assess function, line, and branch coverage from actual executions. Every
reachable uncovered Radiation branch must be assessed for a meaningful fixture
or local test. Record a concrete reason for each remaining gap, including build
configuration and source location. Absence of a matching THC operation alone is
not a reason to omit a reachable local check. No coverage percentage is presumed.
