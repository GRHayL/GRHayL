# Reconstruction Tests And Reference Data

> PLM, WENO-z, and ET_Legacy PPM test binaries, reference-data generators, stencil/dimension guards, and CI versus `run_tests.sh` selection. · Status: confirmed
> Up: [Reconstruction](index.md)

## Summary

Reconstruction tests comprise `unit_test_PLM_reconstruction.c`, `unit_test_WENOZ_reconstruction.c`, and `unit_test_ET_Legacy_reconstruction.c` (which provides GRHayL's only PPM coverage). PLM and WENO-z generators evaluate the same public wrappers as replay, whereas ET_Legacy reconstruction compares against externally supplied legacy outputs and never consumes its locally generated perturbed input.

## Detail

### Test Binaries, Reference Files, And Stencil Guards

Claim status: contested; contradiction: CONTR-0012. See [CONTR-0012](../../contradictions.md#contr-0012).

| Test & generator | Binary files | Valid range, guards & limits |
| --- | --- | --- |
| `unit_test_PLM_reconstruction.c` / `data_gen/unit_test_data_PLM_reconstruction.c` | `PLM_reconstruction_{input,output,output_pert}.bin` | Minmod, MC, superbee over `2 <= index < arraylength - 1` (requires length `>= 4`); boundary slots zeroed by generator and ignored by replay; implementation-derived regression. (Ev: `Unit_Tests/unit_test_PLM_reconstruction.c`) |
| `unit_test_WENOZ_reconstruction.c` / `data_gen/unit_test_data_WENOZ_reconstruction.c` | `WENOZ_reconstruction_{input,output,output_pert}.bin` | `ghl_wenoz_reconstruction` over `3 <= index < arraylength - 2` (requires length `>= 6`); boundary slots zeroed; 5-point helper untested directly; does not prove positivity. (Ev: `Unit_Tests/unit_test_WENOZ_reconstruction.c`) |
| `unit_test_ET_Legacy_reconstruction.c` / `data_gen/unit_test_data_ET_Legacy_reconstruction.c` | `ET_Legacy_reconstruction_{input,output,output_pert}.bin` | 3-direction density/pressure/velocity PPM replay (requires direction length `>= 8` and checked cubic allocation); default PPM params and `Gamma = 2`; generator writes unused `*_input_pert.bin`; no dedicated PPM test/generator exists. (Ev: `Unit_Tests/unit_test_ET_Legacy_reconstruction.c`) |

Each replay binary branches on `ghl_pert_test_fail` and calls `ghl_error` on mismatch, rejecting missing binary headers and undersized stencils.

### Runner Versus Workflow Selection

- `.github/run_tests.sh` downloads `ET_Legacy_reconstruction_*` and `PLM_reconstruction_*` via `.github/et-legacy-testdata-ref` and runs both; it does **not** run `unit_test_WENOZ_reconstruction`.
- Every compiler workflow in `.github/workflows/` runs `PLM_reconstruction` and `WENOZ_reconstruction` in its `reconstruction` matrix (downloading from `TestData` `main`) and runs `reconstruction` in its `ET-Legacy` matrix (using `.github/et-legacy-testdata-ref`).

Evidence (CI behavior): `.github/run_tests.sh::unit_test_PLM_reconstruction`; corroboration: `.github/workflows/github-actions-Ubuntu-gcc.yml::WENOZ_reconstruction`.

## Sources

- [Unit_Tests/unit_test_PLM_reconstruction.c](../../../Unit_Tests/unit_test_PLM_reconstruction.c) - PLM replay test
- [Unit_Tests/data_gen/unit_test_data_PLM_reconstruction.c](../../../Unit_Tests/data_gen/unit_test_data_PLM_reconstruction.c) - PLM generator
- [Unit_Tests/unit_test_WENOZ_reconstruction.c](../../../Unit_Tests/unit_test_WENOZ_reconstruction.c) - WENO-z replay test
- [Unit_Tests/data_gen/unit_test_data_WENOZ_reconstruction.c](../../../Unit_Tests/data_gen/unit_test_data_WENOZ_reconstruction.c) - WENO-z generator
- [Unit_Tests/unit_test_ET_Legacy_reconstruction.c](../../../Unit_Tests/unit_test_ET_Legacy_reconstruction.c) - ET_Legacy PPM replay test
- [Unit_Tests/data_gen/unit_test_data_ET_Legacy_reconstruction.c](../../../Unit_Tests/data_gen/unit_test_data_ET_Legacy_reconstruction.c) - ET_Legacy reconstruction input generator
- [.github/run_tests.sh](../../../.github/run_tests.sh) - local runner selection
- [.github/workflows/github-actions-Ubuntu-gcc.yml](../../../.github/workflows/github-actions-Ubuntu-gcc.yml) - CI `reconstruction` and `ET-Legacy` jobs

## See Also

- Depends on: [Face And Stencil Conventions](face-and-stencil-conventions.md)
- Depends on: [ET Legacy Comparison](../../tests/et-legacy-comparison.md)
- See also: [Runner And Build Outputs](../../tests/runner-and-build-outputs.md)
- See also: [Contradictions](../../contradictions.md#contr-0012)
