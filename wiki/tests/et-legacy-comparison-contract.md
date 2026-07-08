# ET Legacy Comparison Contract

This page owns Unit_Tests routing for `ET_Legacy_*` replay comparisons. These
tests are upstream GRHayL behavior regressions against legacy comparison data.
They are not direct GRHayLib/Cactus verification; implementation-only Cactus
thorn checks route through
[GRHayLib verification and drift](../implementations/grhaylib/verification-and-drift.md).

## Source Authority

- ET Legacy test sources: `Unit_Tests/unit_test_ET_Legacy_*.c`.
- ET Legacy input-side generators: `Unit_Tests/data_gen/unit_test_data_ET_Legacy_*.c`.
- Normal fixture download behavior: [run_tests.sh](../../.github/run_tests.sh)
  and [.github/workflows](../../.github/workflows/).
- Broad docs gap wording, read-only evidence: [mainpage.md](../../docs/raw/mainpage.md).
- Naming/provenance caveats: [contradictions.md](../contradictions.md).

## Normal Run Semantics

Normal scripted/workflow runs download `ET_Legacy_*_input.bin`,
`ET_Legacy_*_output.bin`, and `ET_Legacy_*_output_pert.bin` from TestData before
running each comparison. Repo-visible generators create input-side files such as
`ET_Legacy_*_input.bin` and `ET_Legacy_*_input_pert.bin`; they do not create
trusted output for all ET Legacy families. Do not infer from generator presence
that normal runs regenerate trusted legacy outputs.

The test sources read downloaded trusted output and perturbed output, compute
current GRHayL results, then compare with `ghl_pert_test_fail` or related helper
wrappers. Exact binary ordering stays in paired test/generator source files.

## Test Families

| Test | Behavior area | Fixture family | Generator evidence | Owner route |
| --- | --- | --- | --- | --- |
| `unit_test_ET_Legacy_conservs.c` | Primitive-to-conservative conversion plus stress-energy legacy replay. | `ET_Legacy_conservs_{input,output,output_pert}.bin` | `unit_test_data_ET_Legacy_conservs.c` writes input and perturbed input; source comments say output comes from IllinoisGRMHD/ET legacy path. | [Con2Prim tests and fixtures](../gems/con2prim/tests-and-fixtures.md), [Core tests and fixtures](../core/tests-and-fixtures.md) |
| `unit_test_ET_Legacy_primitives.c` | Conservative-to-primitive recovery legacy replay. | `ET_Legacy_primitives_{input,output,output_pert}.bin` | `unit_test_data_ET_Legacy_primitives.c` writes input and perturbed input; test reads trusted and perturbed output. | [Con2Prim tests and fixtures](../gems/con2prim/tests-and-fixtures.md) |
| `unit_test_ET_Legacy_flux_source.c` | Flux/source RHS legacy replay. | `ET_Legacy_flux_source_{input,output,output_pert}.bin` | `unit_test_data_ET_Legacy_flux_source.c` writes input and perturbed input; test reads trusted and perturbed RHS outputs. | [Flux Source source terms](../gems/flux-source/source-terms-contract.md) |
| `unit_test_ET_Legacy_HLL_flux.c` | Induction vector-potential HLL flux legacy replay. | `ET_Legacy_HLL_flux_{input,output,output_pert}.bin` | `unit_test_data_ET_Legacy_HLL_flux.c` writes input and perturbed input; test reads trusted vector-potential output. | [Induction HLL flux contract](../gems/induction/hll-flux-contract.md), [Induction tests and fixtures](../gems/induction/tests-and-fixtures.md) |
| `unit_test_ET_Legacy_induction_gauge_rhs.c` | BSSN interpolation plus induction vector-potential/gauge RHS replay. | `ET_Legacy_induction_gauge_rhs_{input,output,output_pert}.bin` | `unit_test_data_ET_Legacy_induction_gauge_rhs.c` writes input and perturbed input; test reads trusted RHS output. | [Induction gauge RHS contract](../gems/induction/gauge-rhs-contract.md), [Induction tests and fixtures](../gems/induction/tests-and-fixtures.md) |
| `unit_test_ET_Legacy_reconstruction.c` | PPM reconstruction legacy replay. | `ET_Legacy_reconstruction_{input,output,output_pert}.bin` | `unit_test_data_ET_Legacy_reconstruction.c` writes input and perturbed input; test reads trusted reconstructed face output. | [Reconstruction tests and fixtures](../gems/reconstruction/tests-and-fixtures.md), [PPM flow](../gems/reconstruction/ppm-flow.md) |

## Cross-Gem Boundaries

- `ET_Legacy_conservs` crosses Con2Prim and Core because it checks primitive
  packing, conservative output, and stress-energy output.
- `ET_Legacy_primitives` is Con2Prim recovery evidence; Core helpers appear as
  support code, not owner of recovery behavior.
- `ET_Legacy_flux_source` routes to Flux_Source source-term/RHS behavior.
- `ET_Legacy_HLL_flux` and `ET_Legacy_induction_gauge_rhs` route to Induction;
  do not classify them as Flux_Source HLLE coverage.
- `ET_Legacy_reconstruction` routes to Reconstruction and PPM behavior.

## GRHayLib/Cactus Caveat

ET Legacy tests run GRHayL unit-test binaries. They can support confidence that
upstream GRHayL behavior still matches legacy comparison fixtures, but they do
not verify `implementations/GRHayLib/` build lists, Cactus parameter parsing,
schedule timing, global lifecycle, or Einstein Toolkit runtime behavior.

## Naming And Provenance Caveats

- `TestData` vs `Test_Data` naming drift routes to
  [contradictions.md](../contradictions.md).
- ET Legacy provenance wording about TestData, IllinoisGRMHD, `GRHayLTestPatch`,
  and Einstein Toolkit should stay source-backed. Use repo-visible runner/docs
  evidence unless the external repositories have been explicitly verified.
- Perturbation wording drift also routes to [contradictions.md](../contradictions.md).
