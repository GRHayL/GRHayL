# Source Map

> Aggregate ingest coverage and remaining scope for GRHayL source trees. · Status: confirmed
> Up: [GRHayL Knowledge Base](../AGENTS.md)

## Scope

This table tracks ingest coverage across the multi-file source trees registered in [raw/SOURCES.md](../raw/SOURCES.md), following the authority tiers in [KB Rules](SCHEMA.md#citation-rules) and [Claim And Evidence Rules](SCHEMA.md#claim-and-evidence-rules). Reverse file-to-page dependencies are computed dynamically with `rg -l -F '<path>' wiki/`; external source tiers live in [raw/SOURCES.md](../raw/SOURCES.md).

## Aggregate Rows

| Source or aggregate | Glob pattern | Authority tier | Status | Ingest | Uncovered scope |
| --- | --- | --- | --- | --- | --- |
| `grhayl-core-and-headers` | `GRHayL/make.code.defn`, `GRHayL/GRHayL_Core/**`, `GRHayL/include/**` | primary-code | living | ingested | Individual internal helper macros outside public struct, metric, `u0`, stress-energy, EOS dispatch, and error entry points. |
| `grhayl-atmosphere-con2prim-eos` | `GRHayL/Atmosphere/**`, `GRHayL/Con2Prim/**`, `GRHayL/EOS/**` | primary-code | living | ingested | Minor solver-internal iteration constants and individual interpolator line-level branches beyond the tabulated catalog. |
| `grhayl-flux-induction-reconstruction` | `GRHayL/Flux_Source/**`, `GRHayL/Induction/**`, `GRHayL/Reconstruction/**` | primary-code | living | ingested | Generated NRPy expression temporaries inside `GRHayL/Flux_Source/calculate_*.c` beyond the documented generator boundary. |
| `grhayl-neutrinos-and-compose` | `GRHayL/Neutrinos/**`, `tools/compose/**` | primary-code | living | ingested | Generated `tmp_*` common-subexpression sequences inside `GRHayL/Neutrinos/NRPyLeakage/*.c` and internal regularizer helper constants. |
| `grhaylib-cactus-implementation` | `implementations/GRHayLib/**` | primary-code; primary-doc | living | ingested | External Cactus/CarpetX build and runtime execution outside static CCL/C inspection in this repository. |
| `unit-tests-and-sample-tables` | `Unit_Tests/**` | primary-test | living | ingested | Individual per-row random seeds inside `Unit_Tests/data_gen/*.c` beyond [Data Generators And Helpers](tests/data-generators-and-helpers.md). |
| `doxygen-raw-docs` | `docs/raw/**` | primary-doc | living | ingested | Styling stylesheets (`docs/raw/custom.css`, `docs/raw/doxygen-awesome.css`) and verbatim `@defgroup` Doxygen blocks linked rather than copied. |
| `scripts-and-ci-automation` | `scripts/**`, `.github/**` | primary-code; primary-test | living | ingested | Runner-image OS patch revisions and GitHub Actions runner internals external to `.github/workflows/*.yml`, `.github/actions/**`, and `.github/run_tests.sh`. |
