# Radiation M1 Compatibility Evidence Boundary

This page limits claims to the library operations present in the checkout.
The neutrino M1 implementation has one primary numerical path, with
a flagged Eulerian Minerbo admissibility fallback for finite candidates that
fail pressure or Eulerian covariance positivity, or fail tensor validation at
exact zero flux. The pointwise operations and source-update contract below
describe that path; they do not establish a downstream grid evolution.

Within the retained pointwise corpus, each named record whose THC reference
output fails the runtime pressure-trace or Eulerian covariance positivity
guard is classified by name as a fallback-policy local check, not a cross-code
replay verdict. For
`rngpkt-v2-rd-radiation-anchor-realizability-cone-a01`, whose endpoints share
the identity spatial metric and an x-aligned flux, the test classifies the
retained reference with a negative-Cxx/E witness computed from the reference
independently of the production covariance validator: current GRHayL publishes
its finite Eulerian covariance admissibility fallback, with the analytic
Eulerian Minerbo fallback pressure and closure factor, instead of reproducing
the nonrealizable reference pressure. The species source records derived from
this anchor are classified the same way, because the corrected current closure
fallback changes the comoving moments feeding the frozen-rate source
evaluation; their comparison uses the named operation under the published
closure fallback, with finite sources and the established radiation/matter
exchange check. Records outside
the named fallback-policy classifications retain strict paired cross-code
comparison; no trusted expected value is refreshed from current GRHayL
output.

## What is testable here

- component order `{N,E,Fx,Fy,Fz}`;
- four-point limiter, sawtooth detection, opacity suppression, and high/low
  blending;
- exactly-once face densitization and the canonical no-cap/no-separate-
  diffusion transport boundary;
- frozen-rate source branches and transactional no-update behavior;
- separate total-number and charged-current lepton exchange.

The corresponding implementations are listed in
[`Radiation_traceability.md`](../../../docs/raw/Radiation_traceability.md). The complete
test, fixture, runner, and CI inventory is in
[M1 tests and fixtures](tests-and-fixtures.md). `configure` discovers the
scoped test sources as ordinary `unit_test_*.c` targets, subject to its HDF5
filtering;
The ordinary `.github/run_tests.sh` executes the M1 tests and downloads the
pinned published fixtures by default; executable-level replay is selected
through explicit fixture arguments, and the CI action retries the pinned
download.

Configured execution and a local run remain scoped evidence rather than
downstream host or physical-validation evidence.

## What is not claimed

Repository-local tests can replay inputs against THC_M1 discrete-operation
outputs, as described in [M1 tests and fixtures](tests-and-fixtures.md).
The ordinary runner acquires published `GRHayL/TestData` Radiation members
selected by [radiation-testdata-ref](../../../.github/radiation-testdata-ref),
or consumes a supplied package through `M1_FIXTURE_DIR`. Fixture layout and
producer provenance route through the tests-and-fixtures leaf and
[production record](../../../docs/raw/Radiation_fresh_fixture_production.md).
A missing selected revision or incomplete package fails acquisition; replay
is not silently skipped.
Replay compares the supplied
records; it does not execute THC_M1 live. It does not establish downstream
framework execution, grid evolution,
schedule/AMR behavior, continuum or complete-evolution equivalence, or
physical validation. A downstream consumer owns any such campaign and its
reporting.

## Canonical-method boundary

The host must use the four-point transport operation for neutrino M1 faces and
provide uncapped metric light-cone speeds. The operation does not provide a
separate diffusion correction. The public diffusion helper is tested
separately and is not selected by this canonical route. The existing
finite-difference/Newton source solver remains the local implicit solver.
The public source policy can opt into branched compatibility updates, which
select explicit thin, thick-equilibrium, or scattering-dominated source
algorithms when their conditions hold. The default policy uses the implicit
solver.
The host must apply one limiter scalar to every coupled species, matter, and
`Y_e` increment; no host implementation is included here.
