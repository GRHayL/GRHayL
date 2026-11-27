# ET Legacy Comparison

> Replay rules, input-generator limits, velocity-limiter compatibility bounds, and downstream proof limits for `unit_test_ET_Legacy_*.c`. · Status: contested
> Up: [Tests](index.md)

## Summary

The `Unit_Tests/unit_test_ET_Legacy_*.c` binaries compare standalone GRHayL C routines against downloaded legacy `ET_Legacy_*_{input,output,output_pert}.bin` binaries from `GRHayL/TestData`. Local `unit_test_data_ET_Legacy_*.c` programs generate inputs only (not legacy outputs), and ET_Legacy passes never verify the `implementations/GRHayLib/` Cactus thorn.

## Detail

### Replay Semantics And Special Compatibility Bounds

Claim status: contested; contradiction: CONTR-0003. See [CONTR-0003](../contradictions.md#contr-0003).

Claim status: contested; contradiction: CONTR-0004. See [CONTR-0004](../contradictions.md#contr-0004).

Claim status: contested; contradiction: CONTR-0012. See [CONTR-0012](../contradictions.md#contr-0012).

Scripted and workflow runs download `ET_Legacy_*_{input,output,output_pert}.bin` and compare GRHayL outputs via `ghl_pert_test_fail` or typed wrappers:
- **`ET_Legacy_conservs`:** its input generator targets ADM Lorentz factors in `1.2`–`5` (below `W_max = 10`) and rejects any row where speed limiting fires, keeping the IllinoisGRMHD comparison off the speed-limiter branch while direct Core limiter tests cover the cap boundary.
- **`ET_Legacy_primitives`:** retains unchanged legacy outputs (which omit an independent `eps` field) and uses a fixed velocity-limiter compatibility allowance derived from the older limiter interval for that intentional behavior change rather than rescaling expectations by the recovered Lorentz factor; pressure comparisons stay independent.
- **`ET_Legacy_induction_gauge_rhs`:** enforces strict nonfinite rejection throughout the compared interior against externally supplied legacy outputs.
- **`ET_Legacy_reconstruction`:** exercises PPM reconstruction paths; its local generator writes perturbed input that replay never reads (`CONTR-0012`).

Evidence (descriptive behavior): `Unit_Tests/unit_test_ET_Legacy_primitives.c`; corroboration: `Unit_Tests/data_gen/unit_test_data_ET_Legacy_conservs.c` and `Unit_Tests/unit_test_ET_Legacy_reconstruction.c`.

### ET Legacy Family Map And GRHayLib Proof Boundary

| Test binary | Target subsystem | Defining module test leaf |
| --- | --- | --- |
| `unit_test_ET_Legacy_conservs.c` | Primitive-to-conservative + stress-energy (Con2Prim & Core) | [Con2Prim Tests](../gems/con2prim/tests-and-reference-data.md), [Core Tests](../core/tests-and-reference-data.md) |
| `unit_test_ET_Legacy_primitives.c` | Con2Prim hybrid recovery | [Con2Prim Tests](../gems/con2prim/tests-and-reference-data.md) |
| `unit_test_ET_Legacy_flux_source.c` | Geometric source terms / flux RHS (`Flux_Source`) | [Flux_Source Tests](../gems/flux-source/tests-and-reference-data.md) |
| `unit_test_ET_Legacy_HLL_flux.c` | Induction vector-potential HLL flux (not Flux_Source HLLE) | [Induction Tests](../gems/induction/tests-and-reference-data.md) |
| `unit_test_ET_Legacy_induction_gauge_rhs.c` | BSSN interpolation + vector-potential/gauge RHS | [Induction Tests](../gems/induction/tests-and-reference-data.md) |
| `unit_test_ET_Legacy_reconstruction.c` | PPM legacy reconstruction | [Reconstruction Tests](../gems/reconstruction/tests-and-reference-data.md) |

ET_Legacy binaries execute standalone GRHayL library code only; they do not test `implementations/GRHayLib/` CCL files, Cactus parameter parsing, scheduling, or runtime cleanup ([GRHayLib Verification](../implementations/grhaylib-verification.md)).

## Sources

- [Unit_Tests/unit_test_ET_Legacy_conservs.c](../../Unit_Tests/unit_test_ET_Legacy_conservs.c) - conservative/stress-energy replay
- [Unit_Tests/data_gen/unit_test_data_ET_Legacy_conservs.c](../../Unit_Tests/data_gen/unit_test_data_ET_Legacy_conservs.c) - Lorentz-factor `1.2`–`5` input generator
- [Unit_Tests/unit_test_ET_Legacy_primitives.c](../../Unit_Tests/unit_test_ET_Legacy_primitives.c) - Con2Prim replay and fixed velocity-limiter compatibility bound
- [Unit_Tests/unit_test_ET_Legacy_flux_source.c](../../Unit_Tests/unit_test_ET_Legacy_flux_source.c) - Flux_Source legacy replay
- [Unit_Tests/unit_test_ET_Legacy_HLL_flux.c](../../Unit_Tests/unit_test_ET_Legacy_HLL_flux.c) - Induction HLL legacy replay
- [Unit_Tests/unit_test_ET_Legacy_induction_gauge_rhs.c](../../Unit_Tests/unit_test_ET_Legacy_induction_gauge_rhs.c) - Induction gauge RHS legacy replay
- [Unit_Tests/unit_test_ET_Legacy_reconstruction.c](../../Unit_Tests/unit_test_ET_Legacy_reconstruction.c) - PPM legacy replay
- [docs/raw/mainpage.md](../../docs/raw/mainpage.md) - `GRHayL/TestData` documentation

## See Also

- Depends on: [Reference Data And Test Drivers](reference-data-and-test-drivers.md)
- Depends on: [Runner And Build Outputs](runner-and-build-outputs.md)
- Contrasts with: [GRHayLib Verification](../implementations/grhaylib-verification.md)
- See also: [Contradictions](../contradictions.md#contr-0012)
