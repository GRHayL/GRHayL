# Generated Boundaries

> Separation between hand-authored source files and generated Doxygen, Makefile, test binary, downloaded data, and NRPy-derived C outputs. · Status: confirmed
> Up: [Architecture](index.md)

## Summary

GRHayL distinguishes hand-authored sources (`GRHayL/`, `docs/raw/`, `implementations/GRHayLib/`, `Unit_Tests/`, `tools/compose/`) from generated build and documentation outputs (`Makefile`, `build/`, `test/`, `lib/`, `docs/html/`, `docs/latex/`), downloaded reference data, and NRPy-derived C code in `Flux_Source/` and `Neutrinos/NRPyLeakage/`. Always validate Doxygen in a unique temporary `OUTPUT_DIRECTORY`, verify staged `Flux_Source` generator output before replacing its checked-in C, and edit `NRPyLeakage` C as described in [Generator Provenance](../gems/neutrinos/generator-provenance.md).

## Detail

### Doxygen Input, Layout, And Temporary Validation

`Doxyfile` sets `INPUT = GRHayL docs/raw`, `RECURSIVE = YES`, `EXCLUDE = GRHayL/include/ghl_unit_tests.h GRHayL/Con2Prim/Tabulated/con2prim_CerdaDuran3D.cc`, `EXCLUDE_PATTERNS = *.py`, `USE_MDFILE_AS_MAINPAGE = docs/raw/mainpage.md`, `OUTPUT_DIRECTORY = docs`, HTML enabled, and LaTeX disabled. Root `DoxygenLayout.xml` (header notes Doxygen 1.8.17) is auto-parsed when present; review it alongside `Doxyfile` whenever docs layout changes.

Because `OUTPUT_DIRECTORY = docs` encloses the tracked source subtree `docs/raw/`, never run repository `Doxyfile` in place for validation. Redirect output to a unique temporary directory while keeping repo-root-relative inputs:

```sh
(
  audit_tmp=$(mktemp -d "${TMPDIR:-/tmp}/grhayl-doxygen.XXXXXX") || exit 1
  trap 'rm -rf -- "$audit_tmp"' EXIT
  { cat Doxyfile; printf '\nOUTPUT_DIRECTORY = "%s/output"\n' "$audit_tmp"; } | doxygen -
)
```

Compare Doxygen warnings only across identical Doxygen versions and configurations; record `doxygen --version` in task notes, never in `wiki/`.

Evidence (descriptive behavior): `Doxyfile::OUTPUT_DIRECTORY`; corroboration: `DoxygenLayout.xml`.

### Build, Cleanup, And GRHayLib Symlink Boundaries

| Generator / path | Output paths | Cleanup & safety boundary |
| --- | --- | --- |
| `configure`, then `make` | `configure` writes `Makefile`, `<builddir>/` (default `build/`; `.check`, header symlinks), `test/`, and `test/data_gen/`; `make` adds objects and `<builddir>/lib/libghl*` | `make clean` deletes library objects/libs and normal test binaries/helper objects while keeping `TEXEOBJS` and `data_gen` files; `make realclean` removes `<builddir>/` (default `build/`), `test/`, and `Makefile`; neither removes downloaded root data files. (Ev: `configure::realclean`) |
| `generate_makefile.sh` | `Makefile`, `build/`, `lib/libgrhayl.a` (broken) | Emits GRHayLib source and invalid `-I./include`; reproduce only in a disposable checkout and never run `make`. (Ev: `generate_makefile.sh::libgrhayl.a`) |
| `implementations/GRHayLib/src/` | Tracked `make.code.defn` + symlinks to `GRHayL/` modules and `GRHayL/include` | Repository inputs, not generated copies; see [GRHayLib Cactus Build](../implementations/grhaylib-cactus-build.md). (Ev: `implementations/GRHayLib/src/make.code.defn::SUBDIRS`) |

Claim status: contested; contradiction: CONTR-0002. See [CONTR-0002](../contradictions.md#contr-0002).

### Test Binaries, Downloaded Data, And Sample Tables

`configure` discovers `Unit_Tests/unit_test_*.c` and `Unit_Tests/data_gen/unit_test_data_*.c` and writes `Makefile` rules that `make tests` and `make datagen` build into `test/` and `test/data_gen/`. `.github/run_tests.sh` downloads `*.bin` reference files, `simple_table.h5`, and compressed StellarCollapse tables (`LS220`, `SLy4`) into repository root and cleans tracked paths on `EXIT`; see [Runner And Build Outputs](../tests/runner-and-build-outputs.md). In HDF5 mode, `unit_test_c2p_nn_guess` creates and removes its own private temporary directory for model files. `Unit_Tests/sample_table/generate_simple_table.py` synthesizes `simple_table.h5`, while the checked-in Hempel SFHo simplified `.h5` file is a tracked sample table, not a `configure` output.

Evidence (descriptive behavior): `Unit_Tests/sample_table/generate_simple_table.py`; corroboration: `Unit_Tests/unit_test_c2p_nn_guess.c`.

### NRPy-Derived C Code In `Flux_Source` And `Neutrinos`

- **`Flux_Source`:** `GRHayL/Flux_Source/generate_flux_source.py` is the symbolic NRPy 2 driver (`generate_source`, `generate_speeds`, `generate_fluxes`) producing `ghl_calculate_source_terms.c`, `ghl_calculate_characteristic_speed_dirn{0,1,2}.c`, and each HLLE variant subdirectory (`hybrid`, `hybrid_entropy`, `tabulated`, `tabulated_entropy`). `generate_flux_source.sh <empty-staging-dir>` clones the pinned NRPy 2 revision and runs the driver into an empty directory outside `GRHayL/Flux_Source`, checking its file list against the `make.code.defn` manifests. See [Generated NRPy Boundary](../gems/flux-source/generated-nrpy-boundary.md).
- **`Neutrinos/NRPyLeakage`:** Checked-in C files and headers govern compiled behavior, including hand-applied C corrections not regenerated from ancestral notebooks; see [Generator Provenance](../gems/neutrinos/generator-provenance.md) and [Physics Model And Sources](../gems/neutrinos/physics-model-and-sources.md).

Evidence (generated evidence): `GRHayL/Flux_Source/generate_flux_source.py::generate_fluxes`; corroboration: `GRHayL/Flux_Source/generate_flux_source.sh`.

## Sources

- [Doxyfile](../../Doxyfile) - `INPUT`, `EXCLUDE`, `EXCLUDE_PATTERNS`, `USE_MDFILE_AS_MAINPAGE`, `OUTPUT_DIRECTORY`
- [DoxygenLayout.xml](../../DoxygenLayout.xml) - auto-parsed Doxygen layout file
- [configure](../../configure) - generated build tree, `clean`, and `realclean` targets
- [generate_makefile.sh](../../generate_makefile.sh) - broken legacy static `Makefile` generator
- [implementations/GRHayLib/src/make.code.defn](../../implementations/GRHayLib/src/make.code.defn) - Cactus thorn `SUBDIRS` manifest
- [GRHayL/Flux_Source/generate_flux_source.py](../../GRHayL/Flux_Source/generate_flux_source.py) - `generate_source`, `generate_speeds`, `generate_fluxes`
- [GRHayL/Flux_Source/generate_flux_source.sh](../../GRHayL/Flux_Source/generate_flux_source.sh) - pinned NRPy 2 fetch and staging wrapper
- [Unit_Tests/sample_table/generate_simple_table.py](../../Unit_Tests/sample_table/generate_simple_table.py) - analytic `simple_table.h5` generator
- [Unit_Tests/unit_test_c2p_nn_guess.c](../../Unit_Tests/unit_test_c2p_nn_guess.c) - private temporary directory for HDF5 model tests

## See Also

- Depends on: [Build And Configure](build-and-configure.md)
- Implements: [Generated NRPy Boundary](../gems/flux-source/generated-nrpy-boundary.md)
- See also: [Generator Provenance](../gems/neutrinos/generator-provenance.md)
- See also: [Runner And Build Outputs](../tests/runner-and-build-outputs.md)
