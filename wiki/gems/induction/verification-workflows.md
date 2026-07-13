# Induction Verification Workflows

This page routes build, targeted local runs, and CI evidence for Induction
verification. Fixture ownership and exact fixture families belong in
[tests and fixtures](tests-and-fixtures.md); this page only says which commands
to run after the normal build and fixture placement route is satisfied.

Commands below are supported repository routes, not claims of current
execution. Keep configured, compiled-unrun, locally executed, and CI-configured
statuses separate.

## Build Route

Use the repo-wide configure/build workflow before targeted Induction replay:

```bash
./configure -r
make tests
export LD_LIBRARY_PATH="$(pwd)/build/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
```

The full local CI driver also exports `LD_LIBRARY_PATH` with `build/lib` before
running tests. The GitHub workflow snippets use the installed `lib` directory
or platform-specific setup after their compile action. For repo-wide configure
flags, HDF5 probing, and CI compile details, see [build and CI](../../build-and-ci.md),
[workflows](../../workflows.md), and [`configure`](../../../configure).

All Induction library sources and named tests are build-configured in both
default and `--disable-hdf5` modes. This page does not upgrade manifest or
configured-target evidence into a compile pass.

## Full Local Route

Run the full local-style sequence when an Induction change may interact with
other GRHayL modules:

```bash
.github/run_tests.sh
```

That script runs `./configure -r`, runs `make tests`, sets the runtime library
path, downloads binary fixtures and EOS tables needed by the full suite, runs
all compiled tests it owns, then removes downloaded root-level fixture files.
It ends with `rm -f ./*.bin ./*.h5 ./*.bz2`; run it only in a clean disposable
worktree, never beside caller-owned root-level data with those suffixes.
Use [tests and fixtures](tests-and-fixtures.md) for Induction fixture names
instead of duplicating them here.

## Targeted Local Routes

After normal build and fixture placement, use these smaller checks.

HLL vector-potential flux:

```bash
./test/unit_test_HLL_flux
./test/unit_test_ET_Legacy_HLL_flux
```

Interpolator variants:

```bash
./test/unit_test_induction_ccc_ADM
./test/unit_test_induction_ccc_BSSN
./test/unit_test_induction_vvv_ADM
```

Gauge RHS and ET Legacy compatibility:

```bash
./test/unit_test_ET_Legacy_induction_gauge_rhs
```

Do not regenerate fixtures for normal verification; use the fixture download or
placement route documented by [tests and fixtures](tests-and-fixtures.md) and
[`.github/run_tests.sh`](../../../.github/run_tests.sh).

## CI Mapping

Repo-visible workflow evidence maps Induction coverage this way:

- `ET-Legacy` matrix entries include `induction_gauge_rhs` and `HLL_flux`, so
  they run `unit_test_ET_Legacy_induction_gauge_rhs` and
  `unit_test_ET_Legacy_HLL_flux` with ET Legacy fixture triplets.
- `induction-interpolators` jobs download induction interpolation fixtures and
  run `./test/unit_test_induction_${{ matrix.centering }}_${{ matrix.metric }}`
  for `ccc_ADM`, `ccc_BSSN`, and `vvv_ADM`; visible matrices exclude
  `vvv_BSSN`.
- `induction-flux` jobs download `HLL_flux` fixture files and run
  `./test/unit_test_HLL_flux`.

These workflow files are evidence only; do not edit CI scripts or workflow
definitions for a KB-only verification-route update.

All five workflow files carry these jobs. Push and pull-request triggers ignore
`wiki/**` and Markdown-only changes, so KB-only updates are workflow-only
configuration evidence rather than an automatic Induction run. Scheduled
triggers remain independent of path filters.

## Legacy Generator Boundary

The normal build uses `configure` plus `scripts/parser`; that route parses
Induction manifests correctly. The separate legacy
[`generate_makefile.sh`](../../../generate_makefile.sh) route is **broken** for
current Induction manifests. A clean disposable run exits zero but interprets
the comment `# Primary make.code.defn for GRHayL chalice` as source tokens,
emitting bogus targets such as `GRHayL/Induction/#`, `Primary`, `for`, `GRHayL`,
and `chalice`.

Cause is the trailing continuation on the last `SRCS` line in
[`Interpolators/make.code.defn`](../../../GRHayL/Induction/Interpolators/make.code.defn):
the legacy AWK range continues into the next manifest. Do not run `make` from
that generated file. Maintainer choice remains whether to repair legacy parser,
change manifest formatting, or retire the legacy route; normal `configure`
build evidence is separate.

## HDF5 Note

The Induction tests listed on this page have no visible HDF5-only EOS table
dependency. Repo-wide configure/build behavior may still probe HDF5 by default,
and full-suite routes may download EOS tables for other modules.

## Repo-Local References

- [Test map](../../test-map.md)
- [Workflows](../../workflows.md)
- [Build and CI](../../build-and-ci.md)
- [Induction tests and fixtures](tests-and-fixtures.md)
- [`configure`](../../../configure)
- [`.github/run_tests.sh`](../../../.github/run_tests.sh)
- [Ubuntu GCC workflow](../../../.github/workflows/github-actions-Ubuntu-gcc.yml)
- [Ubuntu clang workflow](../../../.github/workflows/github-actions-Ubuntu-clang.yml)
- [Ubuntu Intel workflow](../../../.github/workflows/github-actions-Ubuntu-intel.yml)
- [macOS GCC workflow](../../../.github/workflows/github-actions-MacOS-gcc.yml)
- [macOS clang workflow](../../../.github/workflows/github-actions-MacOS-clang.yml)
