# Induction Verification Workflows

> Targeted local build/replay commands, CI job mapping, and build-manifest notes for `GRHayL/Induction/`. · Status: confirmed
> Up: [Induction](index.md)

## Summary

All Induction sources and unit tests build in both default and `--disable-hdf5` configurations with no HDF5 EOS table dependency. Use targeted binary invocations after `./configure -r && make tests` (inside a disposable tree) or `.github/run_tests.sh` for cross-module checks; CI runs Induction across `ET-Legacy`, `induction-interpolators`, and `induction-flux`.

## Detail

### Targeted Local Commands And CI Job Mapping

Build per [Build And Configure](../../architecture/build-and-configure.md) (`./configure -r && make tests` and `export LD_LIBRARY_PATH="$(pwd)/build/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"`), place required reference binaries from [Induction Tests And Reference Data](tests-and-reference-data.md), and run only the affected binaries:

```bash
# HLL vector-potential flux:
./test/unit_test_HLL_flux
./test/unit_test_ET_Legacy_HLL_flux
# Interpolator variants:
./test/unit_test_induction_ccc_ADM
./test/unit_test_induction_ccc_BSSN
./test/unit_test_induction_vvv_ADM
# Gauge RHS and ET_Legacy compatibility:
./test/unit_test_ET_Legacy_induction_gauge_rhs
```

- **CI workflows:** Every compiler workflow runs `ET-Legacy` (`induction_gauge_rhs`, `HLL_flux`), `induction-interpolators` (`ccc_ADM`, `ccc_BSSN`, `vvv_ADM`; excludes nonexistent `vvv_BSSN`), and `induction-flux` (`unit_test_HLL_flux`).
- **Build-generator note:** `configure` + `scripts/parser` parse Induction manifests cleanly (`GRHayL/Induction/Interpolators/make.code.defn` terminates its final `SRCS` line without trailing continuation); do not run `make` on broken legacy `generate_makefile.sh` output ([Build And Configure](../../architecture/build-and-configure.md)).

Evidence (CI behavior): `github-actions-Ubuntu-gcc.yml::induction-interpolators`; corroboration: `github-actions-Ubuntu-gcc.yml::induction-flux` and `run_tests.sh::unit_test_HLL_flux`.

## Sources

- [.github/run_tests.sh](../../../.github/run_tests.sh) - local runner Induction commands
- [.github/workflows/github-actions-Ubuntu-gcc.yml](../../../.github/workflows/github-actions-Ubuntu-gcc.yml) - `induction-interpolators`, `induction-flux`, and `ET-Legacy` jobs
- [GRHayL/Induction/Interpolators/make.code.defn](../../../GRHayL/Induction/Interpolators/make.code.defn) - interpolator manifest

## See Also

- Depends on: [Induction Tests And Reference Data](tests-and-reference-data.md)
- Depends on: [Build And Configure](../../architecture/build-and-configure.md)
- See also: [CI Workflows](../../architecture/ci-workflows.md)
