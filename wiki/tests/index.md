# Unit_Tests Hub

This page routes cross-cutting `Unit_Tests/` questions below
[Test Map](../test-map.md). Use it for runner behavior, generated artifacts,
fixture lifecycle, HDF5 sample tables, ET Legacy replay, expected-failure keys,
and coverage/gap routing. It is a router, not source authority.

Source authority stays in repo-local tests, runner/build scripts, and existing
KB routers: [Unit_Tests](../../Unit_Tests/), [configure](../../configure),
[.github/run_tests.sh](../../.github/run_tests.sh),
[.github/workflows](../../.github/workflows/), [Test Map](../test-map.md),
[Build And CI](../build-and-ci.md), and
[Generated Boundaries](../generated-boundaries.md). If this page conflicts with
those files, trust repo-local ground truth and update this router.

## Cross-Cutting Pages

| Page | Use it for |
| --- | --- |
| [Runner and generated artifacts](runner-and-generated-artifacts.md) | `./configure -r`, `make tests`, `make datagen`, `LD_LIBRARY_PATH`, `test/`, `test/data_gen/`, root-level downloaded/generated fixtures, cleanup, HDF5 target filtering, and workflow/runner differences. |
| [Fixture lifecycle and harness contract](fixture-lifecycle-and-harness-contract.md) | How fixture files are generated, downloaded, replayed, perturbed, or intentionally refreshed across Unit_Tests. |
| [HDF5 sample tables](hdf5-sample-tables.md) | Checked-in sample table inputs, generated `simple_table.h5`, downloaded stellar-collapse tables, and no-HDF5 gates. |
| [ET Legacy comparison contract](et-legacy-comparison-contract.md) | `unit_test_ET_Legacy_*` replay semantics, trusted legacy output fixtures, and downstream/compatibility caveats. |
| [Expected-failure and error keys](expected-failure-and-error-keys.md) | `unit_test_code_error` key ranges, expected process failures, HDF5-only skips, and error-path routing. |
| [Unit-test coverage and gap matrix](unit-test-coverage-and-gap-matrix.md) | Behavior-family coverage, weak spots, missing direct tests, installed `ghl_unit_tests.h` caveats, and review impact paths. |

## Existing Per-Gem Fixture Pages

These pages own gem-specific fixture inventories. Link to them instead of
duplicating their tables here:

- [Core tests and fixtures](../core/tests-and-fixtures.md)
- [Con2Prim tests and fixtures](../gems/con2prim/tests-and-fixtures.md)
- [EOS tests and fixtures](../gems/eos/tests-and-fixtures.md)
- [Flux_Source tests and fixtures](../gems/flux-source/tests-and-fixtures.md)
- [Induction tests and fixtures](../gems/induction/tests-and-fixtures.md)
- [Induction verification workflows](../gems/induction/verification-workflows.md)
- [Neutrinos tests and fixtures](../gems/neutrinos/tests-and-fixtures.md)
- [Reconstruction tests and fixtures](../gems/reconstruction/tests-and-fixtures.md)
- [GRHayLib verification and drift](../implementations/grhaylib/verification-and-drift.md)

Flux_Source fixture coverage routes through
[Flux_Source tests and fixtures](../gems/flux-source/tests-and-fixtures.md),
[Test Map](../test-map.md), and [Flux Source](../gems/flux-source.md).
That page also owns the explicit exclusion for Induction vector-potential HLL
tests that should not be counted as Flux_Source HLLE coverage.

## Policy

This hub is a router. Keep source, tests, runner scripts, CI workflows,
`README.md`, `configure`, and Doxygen source as authority; update this page when
those authoritative files change or when routing needs clarification.

Do not add source-tracking hashes, source hashes, `mtime`, stored fingerprints,
unnecessary dates, or separate maintenance logs. Git history is the durable
operation log.

## Repo-Local References

- [AGENTS.md](../../AGENTS.md)
- [wiki/index.md](../index.md)
- [wiki/catalog.md](../catalog.md)
- [wiki/test-map.md](../test-map.md)
- [wiki/build-and-ci.md](../build-and-ci.md)
- [wiki/generated-boundaries.md](../generated-boundaries.md)
- [configure](../../configure)
- [.github/run_tests.sh](../../.github/run_tests.sh)
- [.github/workflows](../../.github/workflows/)
