# Source Map

> Glob patterns, authority tiers, and uncovered scope for GRHayL multi-file source aggregates. · Status: confirmed
> Up: [GRHayL Knowledge Base](../AGENTS.md)

## Scope

This table maps the multi-file source trees registered in [raw/SOURCES.md](../raw/SOURCES.md) to glob patterns, authority tiers, and uncovered scope; its `Glob pattern` column alone defines each aggregate's path membership, following [KB Rules](SCHEMA.md#citation-rules) and [Claim And Evidence Rules](SCHEMA.md#claim-and-evidence-rules). Each aggregate's `Status` and `Ingest` and all external source tiers live only in [raw/SOURCES.md](../raw/SOURCES.md). Reverse file-to-page dependencies are computed dynamically with `rg -l -F '<basename>' wiki/`.

## Aggregate Rows

| Source or aggregate | Glob pattern | Authority tier | Uncovered scope |
| --- | --- | --- | --- |
| `grhayl-core-and-headers` | `GRHayL/make.code.defn`, `GRHayL/GRHayL_Core/**`, `GRHayL/include/**` | primary-code | Individual internal helper macros outside public struct, metric, `u0`, stress-energy, EOS dispatch, and error entry points. |
| `grhayl-atmosphere-con2prim-eos` | `GRHayL/Atmosphere/**`, `GRHayL/Con2Prim/**`, `GRHayL/EOS/**` | primary-code | Minor solver-internal iteration constants and individual interpolator line-level branches beyond the tabulated catalog. |
| `grhayl-flux-induction-reconstruction` | `GRHayL/Flux_Source/**`, `GRHayL/Induction/**`, `GRHayL/Reconstruction/**` | primary-code | Generated NRPy expression temporaries inside `GRHayL/Flux_Source/**/ghl_calculate_*.c` beyond the documented generator boundary. |
| `grhayl-neutrinos-and-compose` | `GRHayL/Neutrinos/**`, `tools/compose/**` | primary-code | Generated `tmp_*` common-subexpression sequences inside `GRHayL/Neutrinos/NRPyLeakage/*.c` and internal regularizer helper constants. |
| `grhaylib-cactus-implementation` | `implementations/GRHayLib/**` | primary-code; primary-doc | External Cactus/CarpetX build and runtime execution outside static CCL/C inspection in this repository. |
| `unit-tests-and-sample-tables` | `Unit_Tests/**` | primary-test | Individual per-row random seeds inside `Unit_Tests/data_gen/*.c` beyond [Data Generators And Helpers](tests/data-generators-and-helpers.md). |
| `doxygen-raw-docs` | `docs/raw/**` | primary-doc | Styling stylesheets (`docs/raw/custom.css`, `docs/raw/doxygen-awesome.css`) and verbatim `@defgroup` Doxygen blocks linked rather than copied. |
| `scripts-and-ci-automation` | `scripts/**`, `.github/**` | primary-code; primary-test | Runner-image OS patch revisions and GitHub Actions runner internals external to `.github/workflows/*.yml`, `.github/actions/**`, and `.github/run_tests.sh`. |
