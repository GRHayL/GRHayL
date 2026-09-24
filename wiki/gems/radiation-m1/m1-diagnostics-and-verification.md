# M1 diagnostics and verification boundary

Diagnostics in the current Radiation M1 implementation describe local
numerical decisions and failures. They are essential for interpreting a
downstream run, but they are not a substitute for a mesh-level verification
campaign.

The public structures and flags are declared in
[`ghl_m1.h`](../../../GRHayL/include/ghl_m1.h). The source paths that update
them are included by
[`Radiation/make.code.defn`](../../../GRHayL/Radiation/make.code.defn) and
[`Radiation/Neutrinos/make.code.defn`](../../../GRHayL/Radiation/Neutrinos/make.code.defn).

## Transport diagnostics

`ghl_m1_four_point_transport_diagnostics` reports, for the five components in
`{N,E,Fx,Fy,Fz}`:

- the componentwise limiter `phi`;
- the componentwise sawtooth decision;
- the common face opacity-suppression factor; and
- the maximum adjacent speed used by the Rusanov low flux.

These fields explain why a face was closer to the centered high flux or the
Rusanov low flux. They do not report a closure solve, a global CFL condition,
or a continuum error estimate. The implementation is
[`ghl_m1_four_point_blended_rusanov.c`](../../../GRHayL/Radiation/ghl_m1_four_point_blended_rusanov.c).

## Closure and diffusion diagnostics

The shared `ghl_m1_diagnostics` record can carry the closure reduced flux,
root residual, root iteration count and status, Eulerian flux factor, Minerbo
factor, realizability-repair flag, `Jthick` validity, and optional diffusion
blend factors. A closure's `four_point_compatibility` field distinguishes the
primary full four-dimensional tensor from the finite admissibility fallback;
it is diagnostic metadata, not a runtime method selector.

The optional diffusion helper reports its blend factor when requested. A
no-op caused by a thin optical depth, invalid thick-limit scalar, or invalid
face velocity is not evidence that the canonical four-point transport was
used; the host must know which operation it called.

## Source and solver diagnostics

`ghl_m1_neutrino_source_diagnostics` records the selected source path, whether
a closure fallback was observed, terminal no-update status, and nested
implicit-solver diagnostics. The latter include:

- Newton iteration count;
- line-search backtracks;
- fallback substep count and whether substepping was used;
- maximum and weighted residual norms; and
- solution-path flags for convergence, projection, backtracking, closure
  fallback, substepping, endpoint acceptance, or terminal failure.

`ghl_m1_neutrino_diagnostics` additionally records provider/rate failures,
source convergence and terminal fallback counts, hard source failures,
number and E/F repair counts, limiter reductions, post-update mean-energy
checks, accumulated repair magnitudes, and accepted rate-product underflow
events. Initialize caller-owned diagnostics with
`ghl_m1_neutrino_diagnostics_initialize` before first use.

Diagnostics do not turn an invalid candidate into a valid candidate. Source
and pair outputs remain transactional: on hard failure or terminal no-update,
the source base and zero exchange packet are retained while the applicable
failure/path counters are updated. After required pointer validation, a hard
pair failure increments `source_failures` exactly once per species from the
caller’s pre-call accumulator; a terminal no-update increments
`source_terminal_fallbacks` exactly once per species instead. Other caller-owned
neutrino diagnostic fields remain unchanged on an aborted pair transaction, and
partial independent or pair-stage source diagnostics are not published. A
successful active pair update retains the independent-stage diagnostics and
adds exactly one `source_converged` count per species for the paired E/F stage.
When all pair fields are inactive, the paired entry point remains equivalent to
the independent single-species route.

## Local test inventory

The scoped M1 runner
[`run_m1_tests.sh`](../../../Unit_Tests/run_m1_tests.sh) builds or runs these
focused executables:

- [`unit_test_m1_closure_fallback.c`](../../../Unit_Tests/unit_test_m1_closure_fallback.c)
- [`unit_test_m1_diffusion_flux.c`](../../../Unit_Tests/unit_test_m1_diffusion_flux.c)
- [`unit_test_m1_error_handling.c`](../../../Unit_Tests/unit_test_m1_error_handling.c)
- [`unit_test_m1_fd_jacobian.c`](../../../Unit_Tests/unit_test_m1_fd_jacobian.c)
- [`unit_test_m1_neutrino_rusanov_flux.c`](../../../Unit_Tests/unit_test_m1_neutrino_rusanov_flux.c)
- [`unit_test_m1_neutrino_seeded_invariants.c`](../../../Unit_Tests/unit_test_m1_neutrino_seeded_invariants.c)
- [`unit_test_m1_neutrino_source_update.c`](../../../Unit_Tests/unit_test_m1_neutrino_source_update.c)
- [`unit_test_m1_rate_provider.c`](../../../Unit_Tests/unit_test_m1_rate_provider.c)
- [`unit_test_m1_thcm1_blended_rusanov.c`](../../../Unit_Tests/unit_test_m1_thcm1_blended_rusanov.c)
- [`unit_test_rusanov_flux.c`](../../../Unit_Tests/unit_test_rusanov_flux.c)

The [M1 test guide](../../../Unit_Tests/README.m1.md) documents the scoped
build, generated-provider mode, stored fixtures, and current reference
families. The variable-volume transport fixtures include the prepared-face
operation; the seeded invariant tests include exact-zero and underflow/overflow
state cases; source tests include endpoint and transactional checks.

## What the evidence establishes

These tests provide local evidence for closure/realizability, physical and
Rusanov transport, four-point blending, optional diffusion, rate-provider
validation, finite-difference Jacobians, source endpoints, and error paths.
Stored pointwise and prepared-face comparisons validate the corresponding
discrete operations and their input/output conventions.

They do not by themselves establish:

- a complete finite-volume mesh evolution or stable global timestep policy;
- AMR, boundaries, refluxing, reconstruction, or framework integration;
- continuum convergence or physical validation against a production transport
  code;
- a successful remote CI run merely because a workflow selects these tests; or
- line/branch coverage without inspecting actual execution artifacts.

The stored-reference boundary is maintained in
[`README.m1.md`](../../../Unit_Tests/README.m1.md) and the fixture-family
documents under `Unit_Tests/data/m1_thcm1/`. The package is validated by
`audit_package.py`, and current reference values are not regenerated from the
GRHayL output during a test run.

## Verification practice for downstream hosts

When integrating a host, retain the local diagnostics and record the selected
transport/source path, metric and volume convention, rate-provider status,
repair budgets, limiter reductions, endpoint failures, and terminal no-update
events. Add host-level tests for stencil assembly, volume preparation, flux
divergence, species pairing, common limiting, matter recovery, and stage
ordering. Those checks are outside the shipped GRHayL unit-test boundary and
must not be implied by a passing local executable.
