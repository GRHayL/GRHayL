# Radiation M1 Compatibility Evidence Boundary

This page limits claims to the library operations present in the checkout.
The neutrino M1 implementation has one primary numerical path, with
a flagged Eulerian Minerbo admissibility fallback for finite non-PSD or
invalid exact-zero-flux closure candidates. The pointwise operations and source-update contract below describe
that path; they do not establish a downstream grid evolution.

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
outputs, as described in [M1 tests and fixtures](tests-and-fixtures.md). A
fresh candidate now has new producer records and consumed inputs for every
family and is published as TestData members for default replay. The imported
[`radiation-testdata-ref`](../../../.github/radiation-testdata-ref)
pin selects the TestData revision for the ordinary download route; that
revision must be published before either the runner or the action can acquire
it.
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
