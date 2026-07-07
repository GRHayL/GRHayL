# GRHayL Knowledge Base

This repository carries GRHayL source, tests, implementations, and Doxygen
source documentation. Start with the KB router in `wiki/index.md`, then use
`wiki/catalog.md` for term lookup. Keep direct source, Doxygen, test, and CI
links in view because repo-local files remain the ground truth.

## Router

| Go to | Use it for |
| --- | --- |
| [KB Index](wiki/index.md) | Root router for KB pages, source authority, page contracts, and source/date policy. |
| [KB Catalog](wiki/catalog.md) | Query-routing table for aliases, gems, source paths, docs, and tests. |
| [Source Map](wiki/source-map.md) | Source tree ownership and dependency routing. |
| [Workflows](wiki/workflows.md) | Common build, test, docs, and review workflows for agents. |
| [Test Map](wiki/test-map.md) | Unit test, reference-data, fixture, and CI test routing. |
| [README](README.md) | Project purpose, gem overview, build/install paths, HDF5 configuration, CI testing, and implementation overview. |
| [Doxygen Main Page](docs/raw/mainpage.md) | Documentation source entry point used by `Doxyfile`; links gem groups, derivation, build/install, CI, and implementation notes. |
| [Core API Docs](docs/raw/GRHayL_Core.dox) | Core Doxygen groups for EOS initialization, struct packing/unpacking, metric helpers, and stress-energy helpers. |
| [Gems](GRHayL/) | Source tree for Atmosphere, Con2Prim, EOS, Flux_Source, Induction, Neutrinos, and Reconstruction modules. |
| [Equations](docs/raw/derivation.md) | GRMHD variables, stress-energy tensor, conservative equations, flux/source terms, and magnetic-field rescaling. |
| [Implementations](implementations/GRHayLib/) | GRHayLib downstream integration files. |
| [Testing](Unit_Tests/) | Unit tests, generated reference-data programs, perturbation checks, and sample EOS table. |
| [CI Workflows](.github/workflows/) | GitHub Actions build and test coverage. |
| [Doxygen Configuration](Doxyfile) | Documentation input/output configuration and generated-output boundary. |

## Where Do I Start?

| Task | Read first |
| --- | --- |
| Find a term, alias, module, or owned page | [KB Catalog](wiki/catalog.md), [KB Index](wiki/index.md) |
| Build or install GRHayL | [README](README.md), [configure](configure) |
| Understand GRHayL's modular design | [README](README.md), [Doxygen Main Page](docs/raw/mainpage.md) |
| Understand generated outputs and build products | [Doxyfile](Doxyfile), [docs/raw](docs/raw/) |
| Check style and contribution rules | [.clang-format](.clang-format), [.editorconfig](.editorconfig) |
| Find structs, EOS setup, metric helpers, or stress-energy helpers | [Core API Docs](docs/raw/GRHayL_Core.dox), [public headers](GRHayL/include/) |
| Work on atmosphere prescriptions | [Atmosphere docs](docs/raw/Atmosphere.dox), [Atmosphere source](GRHayL/Atmosphere/) |
| Work on conservative-to-primitive solvers | [Con2Prim docs](docs/raw/Con2Prim.dox), [Con2Prim source](GRHayL/Con2Prim/) |
| Work on equation-of-state routines | [EOS docs](docs/raw/EOS.dox), [EOS source](GRHayL/EOS/) |
| Work on fluxes, source terms, or characteristic speeds | [Flux Source docs](docs/raw/Flux_Source.dox), [Flux Source source](GRHayL/Flux_Source/) |
| Work on vector-potential induction routines | [Induction docs](docs/raw/Induction.dox), [Induction source](GRHayL/Induction/) |
| Work on neutrino leakage support | [Neutrinos source](GRHayL/Neutrinos/), [radiation header](GRHayL/include/ghl_radiation.h) |
| Work on shock-capturing reconstruction | [Reconstruction docs](docs/raw/Reconstruction.dox), [Reconstruction source](GRHayL/Reconstruction/) |
| Understand the GRMHD derivation | [GRMHD Derivation](docs/raw/derivation.md) |
| Connect GRHayL to downstream infrastructure | [GRHayLib implementation](implementations/GRHayLib/) |
| Run unit tests or inspect CI coverage | [Test Map](wiki/test-map.md), [Unit tests](Unit_Tests/), [CI workflows](.github/workflows/), [run_tests.sh](.github/run_tests.sh) |
| Update documentation pages | [docs/raw](docs/raw/), [Doxyfile](Doxyfile) |

## Source-Tracking And Date Policy

These rules bind `AGENTS.md`, `wiki/`, and any future GRHayL KB manifest or
governance pages added to this repository:

- No source-tracking hashes or hashing of sources.
- No `mtime`.
- Avoid KB dates unless they are absolutely necessary. When retained dates are
  necessary, use `MM-DD-YYYY`.
- Do not output KB maintenance notes to a separate maintenance log. This
  repository already lives in git, so commit history records durable operations.

Source drift is handled by dependency-aware review of changed repository paths
and affected documentation pages, not by stored fingerprints.
