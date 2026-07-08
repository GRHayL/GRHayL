# GRHayL Agent KB

This knowledge base routes agents to GRHayL source, tests, implementations,
and Doxygen source documentation without replacing them. Repo-local files are
the authority: prefer `README.md`, `Doxyfile`, `configure`, `docs/raw/`,
`GRHayL/`, `Unit_Tests/`, `.github/workflows/`, and
`implementations/GRHayLib/` before broad searches.

KB pages summarize navigation and review impact. If a KB page conflicts with
source, headers, tests, CI, or Doxygen source, trust the underlying repo files
and update the KB page.

## Router

| Go to | Use it for | Primary ground truth |
| --- | --- | --- |
| [Catalog](catalog.md) | Query routing by term, alias, gem, and workflow. | `README.md`, `docs/raw/`, `GRHayL/`, `Unit_Tests/` |
| [Source Map](source-map.md) | Source tree ownership and dependency routing. | `GRHayL/`, `GRHayL/include/`, `GRHayL/make.code.defn` |
| [Public API Map](public-api-map.md) | Public headers, structs, function families, and callable surface. | `GRHayL/include/`, `docs/raw/GRHayL_Core.dox` |
| [Core/Chalice](core/index.md) | Shared Core structs, pack/unpack helpers, metrics, `u0`, stress-energy, EOS dispatch, errors, IO, debug utilities, and Core fixtures. | `GRHayL/GRHayL_Core/`, `GRHayL/include/ghl.h`, `docs/raw/GRHayL_Core.dox`, `Unit_Tests/unit_test_grhayl_core_test_suite.c` |
| [Gems](gems/index.md) | Gem-by-gem module summaries. | `GRHayL/Atmosphere/`, `GRHayL/Con2Prim/`, `GRHayL/EOS/`, `GRHayL/Flux_Source/`, `GRHayL/Induction/`, `GRHayL/Neutrinos/`, `GRHayL/Reconstruction/` |
| [Build And CI](build-and-ci.md) | Configure, HDF5, Makefile, install, and CI workflow routes. | `README.md`, `configure`, `.github/workflows/`, `.github/run_tests.sh` |
| [Test Map](test-map.md) | Unit tests, generated reference-data programs, sample EOS table, and fixture routing. | `Unit_Tests/`, `.github/run_tests.sh` |
| [Generated Boundaries](generated-boundaries.md) | Generated outputs, Doxygen output, NRPy-derived files, and build products. | `Doxyfile`, `docs/raw/`, `GRHayL/Flux_Source/nrpy/`, `GRHayL/*/make.code.defn` |
| [Change Impact](change-impact.md) | Which docs/tests to inspect after source changes. | `GRHayL/`, `GRHayL/include/`, `docs/raw/`, `Unit_Tests/` |
| [Contradictions](contradictions.md) | Known mismatches and review notes that need source-backed resolution. | Repo-local files named by each contradiction |
| [Workflows](workflows.md) | Common agent workflows for docs, tests, builds, and review. | `README.md`, `configure`, `.github/run_tests.sh` |
| [Physics Variables](physics/variables-and-conventions.md) | Primitive/conservative variables, stress-energy notation, magnetic rescaling. | `docs/raw/derivation.md`, `GRHayL/include/ghl.h` |
| [Evolution Equation Map](physics/evolution-equation-map.md) | Flux/source/evolution equation routes and implementation entry points. | `docs/raw/derivation.md`, `docs/raw/Flux_Source.dox`, `GRHayL/Flux_Source/` |
| [KB Checks](lint/CHECKS.md) | Manual checks for links, policy, page contracts, source authority, and duplication. | `AGENTS.md`, `tasks*.md`, `docs/raw/`, KB pages |

## Where Do I Start?

| Task | Read first |
| --- | --- |
| Find a term or alias | [Catalog](catalog.md) |
| Build or install GRHayL | [Build And CI](build-and-ci.md), `README.md`, `configure` |
| Understand modular design | `README.md`, `docs/raw/mainpage.md`, [Gems](gems/index.md) |
| Work on Core/chalice shared structs, packing, metrics, `u0`, stress-energy, EOS dispatch, errors, IO, or debug utilities | [Core/Chalice](core/index.md), `docs/raw/GRHayL_Core.dox`, `GRHayL/include/ghl.h` |
| Inspect public structs or API entry points | [Public API Map](public-api-map.md), `GRHayL/include/`, `docs/raw/GRHayL_Core.dox` |
| Work on atmosphere prescriptions | [Gems](gems/index.md), `docs/raw/Atmosphere.dox`, `GRHayL/Atmosphere/` |
| Work on conservative-to-primitive solvers | [Gems](gems/index.md), `docs/raw/Con2Prim.dox`, `GRHayL/Con2Prim/` |
| Work on EOS routines or HDF5 behavior | [Gems](gems/index.md), `docs/raw/EOS.dox`, `GRHayL/EOS/`, `configure` |
| Work on fluxes, characteristic speeds, or source terms | [Evolution Equation Map](physics/evolution-equation-map.md), `docs/raw/Flux_Source.dox`, `GRHayL/Flux_Source/` |
| Work on vector-potential induction | [Evolution Equation Map](physics/evolution-equation-map.md), `docs/raw/Induction.dox`, `GRHayL/Induction/` |
| Work on neutrino leakage | [Gems](gems/index.md), `GRHayL/Neutrinos/`, `GRHayL/include/ghl_radiation.h` |
| Work on reconstruction | [Gems](gems/index.md), `docs/raw/Reconstruction.dox`, `GRHayL/Reconstruction/` |
| Understand GRMHD equations | [Physics Variables](physics/variables-and-conventions.md), [Evolution Equation Map](physics/evolution-equation-map.md), `docs/raw/derivation.md` |
| Connect to downstream infrastructure | `implementations/GRHayLib/` |
| Run tests or inspect reference data | [Test Map](test-map.md), `Unit_Tests/`, `.github/run_tests.sh` |
| Check generated-output boundaries | [Generated Boundaries](generated-boundaries.md), `Doxyfile`, `GRHayL/*/make.code.defn` |
| Review likely impact of a change | [Change Impact](change-impact.md), [Contradictions](contradictions.md) |
| Maintain the KB | [KB Checks](lint/CHECKS.md), `AGENTS.md` |

## Source And Date Policy

- No source-tracking hashes or hashing of sources.
- No `mtime` tracking.
- Avoid KB dates unless absolutely necessary. If retained, use `MM-DD-YYYY`.
- Do not write KB maintenance notes to a separate maintenance log.
- Handle source drift by dependency-aware review of changed paths and affected
  pages, not by stored fingerprints.

## Page Contract

Each KB page should:

- State its routing purpose near the top.
- Link to repo-relative ground-truth files for every substantive claim.
- Prefer pointers and concise synthesis over copied Doxygen or source content.
- Include a `Ground Truth References` section only when external web sources
  were used, with official full URLs.
- Avoid source-tracking hashes, `mtime`, stored fingerprints, separate
  maintenance logs, and unnecessary dates.
- Keep links repo-relative and compatible with parallel pages that may be
  created by other agents.
