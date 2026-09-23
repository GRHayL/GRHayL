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
[`TRACEABILITY.md`](../../../GRHayL/Radiation/TRACEABILITY.md). The complete
test, fixture, runner, and CI inventory is in
[M1 tests and fixtures](tests-and-fixtures.md). `configure` discovers the ten
scoped test sources as ordinary `unit_test_*.c` targets, subject to its HDF5
filtering;
`Unit_Tests/run_m1_tests.sh` selects and runs them, and the dedicated Radiation
action invokes that route. The normal `.github/run_tests.sh` path does not
select the scoped M1 runner. The checked-in fixture package is validated by
`Unit_Tests/data/m1_thcm1/audit_package.py` before replay.
Configured execution and a local run remain scoped evidence rather than
downstream host or physical-validation evidence.

## What is not claimed

No external-source comparison, downstream framework run, grid evolution,
schedule/AMR result, complete evolution equivalence, or physical-validation
result is established by the listed library operations. A downstream consumer owns
any such campaign and its reporting.

## Canonical-method boundary

The host must use the four-point transport operation for neutrino M1 faces and
provide uncapped metric light-cone speeds. The operation does not provide a
separate diffusion correction. The public diffusion helper is tested
separately and is not selected by this canonical route. The existing
finite-difference/Newton source solver remains the local implicit solver.
Legacy policy fields may remain for
source or ABI compatibility, but cannot select a different numerical method.
The host must apply one limiter scalar to every coupled species, matter, and
`Y_e` increment; no host implementation is included here.
