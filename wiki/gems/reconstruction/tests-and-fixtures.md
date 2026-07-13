# Reconstruction Tests And Fixtures

## Purpose

This page routes Reconstruction test binaries, data generators, fixture names,
and CI behavior. Use it when checking PLM, WENOZ, ET Legacy reconstruction, or
PPM test coverage. Repo-local tests, generators, `configure`, CI scripts, and
[Test Map](../../test-map.md) remain the authority.

Read with the [Reconstruction gem hub](../reconstruction.md).

## Direct Reconstruction Fixtures

- [Unit_Tests/unit_test_PLM_reconstruction.c](../../../Unit_Tests/unit_test_PLM_reconstruction.c)
  reads `PLM_reconstruction_input.bin`,
  `PLM_reconstruction_output.bin`, and
  `PLM_reconstruction_output_pert.bin`.
  [Unit_Tests/data_gen/unit_test_data_PLM_reconstruction.c](../../../Unit_Tests/data_gen/unit_test_data_PLM_reconstruction.c)
  writes the same fixture names for minmod, MC, and superbee PLM output bars.
- [Unit_Tests/unit_test_WENOZ_reconstruction.c](../../../Unit_Tests/unit_test_WENOZ_reconstruction.c)
  reads `WENOZ_reconstruction_input.bin`,
  `WENOZ_reconstruction_output.bin`, and
  `WENOZ_reconstruction_output_pert.bin`.
  [Unit_Tests/data_gen/unit_test_data_WENOZ_reconstruction.c](../../../Unit_Tests/data_gen/unit_test_data_WENOZ_reconstruction.c)
  writes the same fixture names for WENOZ output bars.
- [Unit_Tests/unit_test_ET_Legacy_reconstruction.c](../../../Unit_Tests/unit_test_ET_Legacy_reconstruction.c)
  reads `ET_Legacy_reconstruction_input.bin`,
  `ET_Legacy_reconstruction_output.bin`, and
  `ET_Legacy_reconstruction_output_pert.bin`.
  [Unit_Tests/data_gen/unit_test_data_ET_Legacy_reconstruction.c](../../../Unit_Tests/data_gen/unit_test_data_ET_Legacy_reconstruction.c)
  writes `ET_Legacy_reconstruction_input.bin` and
  `ET_Legacy_reconstruction_input_pert.bin`; it does not produce either output
  fixture. Replay and CI do not consume the generated perturbed-input file.
  Trusted and perturbed outputs come from downloaded fixture data.

## Method Coverage

- PLM direct coverage lives in `unit_test_PLM_reconstruction.c`; the generator
  produces input, trusted output, and perturbed-output bars for
  `PLM_reconstruction`.
- WENOZ direct coverage lives in `unit_test_WENOZ_reconstruction.c`; the
  generator produces input, trusted output, and perturbed-output bars for
  `WENOZ_reconstruction`.
- ET Legacy reconstruction coverage compares density, pressure, and velocity
  left/right face values across three flux directions using
  `ET_Legacy_reconstruction` fixtures.
- ET Legacy exercises PPM paths; no dedicated `unit_test_PPM_reconstruction.c` or PPM data generator is visible.

All three replay tests use `ghl_pert_test_fail` in an `if` condition and call
`ghl_error` on mismatch. Numerical comparisons can therefore fail these test
binaries. Coverage is still bounded:

- PLM replay checks valid interior stencils, but its generator evaluates every
  array index and accesses outside the allocation at both ends.
- WENOZ replay checks valid interior stencils, but its generator writes full
  output arrays whose boundary entries were never initialized.
- PPM helper routines have only transitive ET Legacy coverage; no dedicated
  PPM or helper test exists.

## Configure Versus Runtime Data

[configure](../../../configure) discovers test and data-generator targets from
`Unit_Tests/unit_test_*.c` and `Unit_Tests/data_gen/unit_test_data_*.c`, then
emits `tests` and `datagen` Makefile targets. That is compile-target
discovery.

Fixture download and test execution are separate runtime behavior. The shell
runner and GitHub Actions jobs download fixture files, place them in the
working directory, and run the compiled test executables.

## CI Routes

`.github/run_tests.sh` runs ET Legacy reconstruction and PLM; workflow reconstruction jobs run PLM and WENOZ; ET-Legacy workflow matrix includes reconstruction.

- [.github/run_tests.sh](../../../.github/run_tests.sh) downloads
  `ET_Legacy_reconstruction_{input,output,output_pert}.bin`, runs
  `./test/unit_test_ET_Legacy_reconstruction`, downloads
  `PLM_reconstruction_{input,output,output_pert}.bin`, then runs
  `./test/unit_test_PLM_reconstruction`.
- [.github/workflows/github-actions-Ubuntu-gcc.yml](../../../.github/workflows/github-actions-Ubuntu-gcc.yml),
  [.github/workflows/github-actions-Ubuntu-clang.yml](../../../.github/workflows/github-actions-Ubuntu-clang.yml),
  [.github/workflows/github-actions-Ubuntu-intel.yml](../../../.github/workflows/github-actions-Ubuntu-intel.yml),
  [.github/workflows/github-actions-MacOS-gcc.yml](../../../.github/workflows/github-actions-MacOS-gcc.yml),
  and
  [.github/workflows/github-actions-MacOS-clang.yml](../../../.github/workflows/github-actions-MacOS-clang.yml)
  define a `reconstruction` job matrix for `PLM_reconstruction` and
  `WENOZ_reconstruction`, and an `ET-Legacy` matrix that includes
  `reconstruction`.

Thus WENOZ is workflow-only relative to the aggregate local runner, while PPM
is reached through ET Legacy in both routes. `run_tests.sh` removes downloaded
`*.bin`, `*.h5`, and `*.bz2` artifacts at its end; individual workflow jobs use
fresh workspaces and contain no corresponding cleanup step.

## Repo-Local References

- [Unit_Tests/unit_test_PLM_reconstruction.c](../../../Unit_Tests/unit_test_PLM_reconstruction.c)
- [Unit_Tests/unit_test_WENOZ_reconstruction.c](../../../Unit_Tests/unit_test_WENOZ_reconstruction.c)
- [Unit_Tests/unit_test_ET_Legacy_reconstruction.c](../../../Unit_Tests/unit_test_ET_Legacy_reconstruction.c)
- [Unit_Tests/data_gen/unit_test_data_PLM_reconstruction.c](../../../Unit_Tests/data_gen/unit_test_data_PLM_reconstruction.c)
- [Unit_Tests/data_gen/unit_test_data_WENOZ_reconstruction.c](../../../Unit_Tests/data_gen/unit_test_data_WENOZ_reconstruction.c)
- [Unit_Tests/data_gen/unit_test_data_ET_Legacy_reconstruction.c](../../../Unit_Tests/data_gen/unit_test_data_ET_Legacy_reconstruction.c)
- [.github/run_tests.sh](../../../.github/run_tests.sh)
- [.github/workflows/](../../../.github/workflows/)
- [configure](../../../configure)
- [wiki/test-map.md](../../test-map.md)
- [wiki/gems/reconstruction.md](../reconstruction.md)
