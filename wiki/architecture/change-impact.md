# Change Impact

> Subsystem map linking edited GRHayL paths to affected documentation, downstream boundaries, and targeted validation checks. · Status: confirmed
> Up: [Architecture](index.md)

## Summary

Use this page after editing source, header, documentation, test, or CI files to select affected documentation routes and targeted verification checks without running unrelated suites. Execute mutating builds or tests inside disposable checkouts, and obtain explicit user permission before adding or expanding tests, assertions, generators, or CI test executions.

## Detail

### Core, Headers, And Build/CI Paths

| Changed path | Affected KB & Doxygen routes | Targeted verification checks |
| --- | --- | --- |
| `GRHayL/GRHayL_Core/` | [Core](../core/index.md), `docs/raw/GRHayL_Core.dox`, `docs/raw/mainpage.md` | `unit_test_grhayl_core_test_suite`, `unit_test_compute_conservs_and_Tmunu`, `unit_test_code_error`, broad `make tests`. (Ev: `GRHayL/GRHayL_Core/make.code.defn`) |
| `GRHayL/include/` | [Public API Map](public-api-map.md), relevant module `.dox`, `README.md`, [GRHayLib Verification](../implementations/grhaylib-verification.md) | `./configure -r && make tests datagen`; check `GRHayLib.h` and `install_headers` parity. (Ev: `GRHayL/include/make.code.defn`) |
| `configure` / `generate_makefile.sh` | [Build And Configure](build-and-configure.md), [Generated Boundaries](generated-boundaries.md), `README.md` | Test default and `--disable-hdf5` `./configure -r && make tests`; inspect `generate_makefile.sh` only in a disposable tree (never run `make`). (Ev: `configure::--disable-hdf5`) |
| `docs/raw/` / `Doxyfile` | [Generated Boundaries](generated-boundaries.md), [CI Workflows](ci-workflows.md) | Run temporary-`OUTPUT_DIRECTORY` Doxygen validation; never generate into tracked `docs/`. (Ev: `Doxyfile::OUTPUT_DIRECTORY`) |
| `.github/` (`workflows/`, `actions/`, `run_tests.sh`) | [CI Workflows](ci-workflows.md), [Runner And Build Outputs](../tests/runner-and-build-outputs.md), [Contradictions](../contradictions.md) | Parse YAML with `BaseLoader`, `bash -n .github/run_tests.sh`, reproduce matching `configure`/`make` locally in disposable tree. (Ev: `run_tests.sh`) |

### Gem Subsystems

| Changed path | Affected KB & Doxygen routes | Targeted verification checks |
| --- | --- | --- |
| `GRHayL/Atmosphere/` | [Atmosphere](../gems/atmosphere.md), `docs/raw/Atmosphere.dox` | `unit_test_grhayl_core_test_suite`, plus indirect Con2Prim/NN atmosphere fallback tests. (Ev: `GRHayL/Atmosphere/make.code.defn`) |
| `GRHayL/Con2Prim/` | [Con2Prim](../gems/con2prim/index.md), `docs/raw/Con2Prim.dox`, [GRHayLib Runtime Parameters](../implementations/grhaylib-runtime-parameters.md) | `unit_test_apply_conservative_limits`, `unit_test_con2prim_multi_method_hybrid`, `unit_test_enforce_primitive_limits_and_compute_u0`, `unit_test_compute_conservs_and_Tmunu`, `unit_test_hybrid_failure`, `unit_test_c2p_nn_guess`, `unit_test_con2prim_tabulated`, ET_Legacy C2P tests. (Ev: `GRHayL/Con2Prim/make.code.defn`) |
| `GRHayL/EOS/` | [EOS](../gems/eos/index.md), `docs/raw/EOS.dox`, [CompOSE Callback](../gems/neutrinos/compose-eos-leakage-callback.md) | `unit_test_piecewise_polytrope`, `unit_test_tabulated_eos`, `unit_test_tabulated_eos_compose`, `unit_test_con2prim_tabulated`, NRPyLeakage tests. Note: `unit_test_tabulated_flux` uses synthetic `test_compute_h_and_cs2.c`, not a table sound-speed oracle. (Ev: `GRHayL/EOS/make.code.defn`) |
| `GRHayL/Flux_Source/` | [Flux_Source](../gems/flux-source/index.md), `docs/raw/Flux_Source.dox`, `docs/raw/derivation.md`, [Generated Boundaries](generated-boundaries.md) | `unit_test_hybrid_flux`, `unit_test_tabulated_flux`, `unit_test_ET_Legacy_flux_source` (HDF5 enabled for tabulated variants). (Ev: `GRHayL/Flux_Source/make.code.defn`) |
| `GRHayL/Induction/` | [Induction](../gems/induction/index.md), `docs/raw/Induction.dox`, `docs/raw/derivation.md` | `unit_test_HLL_flux`, `unit_test_ET_Legacy_HLL_flux`, `unit_test_induction_{ccc_ADM,ccc_BSSN,vvv_ADM}`, `unit_test_ET_Legacy_induction_gauge_rhs`. (Ev: `GRHayL/Induction/make.code.defn`) |
| `GRHayL/Neutrinos/` | [Neutrinos](../gems/neutrinos/index.md), `GRHayL/include/ghl_radiation.h`, [Generated Boundaries](generated-boundaries.md) | `unit_test_nrpyleakage_{physics,classifier_fallback,optically_thin_gas,constant_density_sphere,luminosities}`, `unit_test_code_error`. (Ev: `GRHayL/Neutrinos/make.code.defn`) |
| `GRHayL/Reconstruction/` | [Reconstruction](../gems/reconstruction/index.md); keep `docs/raw/Reconstruction.dox` read-only | `unit_test_PLM_reconstruction`, `unit_test_WENOZ_reconstruction`, `unit_test_ET_Legacy_reconstruction` (PPM tested via ET_Legacy only). (Ev: `GRHayL/Reconstruction/make.code.defn`) |

### Tests, CompOSE Tooling, And GRHayLib Implementation

| Changed path | Affected KB & Doxygen routes | Targeted verification checks |
| --- | --- | --- |
| `Unit_Tests/` & `Unit_Tests/sample_table/` | [Tests](../tests/index.md), [Generated Boundaries](generated-boundaries.md), [Contradictions](../contradictions.md) | Compile (`make tests datagen`) and run changed tests/generators separately with verified data; run `.github/run_tests.sh` only in a disposable checkout. (Ev: `generate_simple_table.py`) |
| `tools/compose/` | [CompOSE Conversion](../gems/neutrinos/compose-eos-conversion-and-regularization.md), [CompOSE Callback](../gems/neutrinos/compose-eos-leakage-callback.md), `tools/compose/README.md` | Run existing Python line/branch coverage and `unit_test_tabulated_eos_compose` C integration checks in disposable output paths. (Ev: `compose_to_grhayl.py`) |
| `implementations/GRHayLib/` | [Implementations](../implementations/index.md), `README.md` | Check CCL, `GRHayLib.h`, symlink, `SUBDIRS`, and parameter/parser parity statically; Cactus build/runtime proof requires an external Cactus/ET environment. (Ev: `implementations/GRHayLib/src/make.code.defn`) |

## Sources

- [GRHayL/GRHayL_Core/make.code.defn](../../GRHayL/GRHayL_Core/make.code.defn) - Core manifest
- [GRHayL/include/make.code.defn](../../GRHayL/include/make.code.defn) - public install-header manifest
- [GRHayL/Atmosphere/make.code.defn](../../GRHayL/Atmosphere/make.code.defn) - Atmosphere manifest
- [GRHayL/Con2Prim/make.code.defn](../../GRHayL/Con2Prim/make.code.defn) - Con2Prim manifest
- [GRHayL/EOS/make.code.defn](../../GRHayL/EOS/make.code.defn) - EOS manifest
- [GRHayL/Flux_Source/make.code.defn](../../GRHayL/Flux_Source/make.code.defn) - Flux_Source manifest
- [GRHayL/Induction/make.code.defn](../../GRHayL/Induction/make.code.defn) - Induction manifest
- [GRHayL/Neutrinos/make.code.defn](../../GRHayL/Neutrinos/make.code.defn) - Neutrinos manifest
- [GRHayL/Reconstruction/make.code.defn](../../GRHayL/Reconstruction/make.code.defn) - Reconstruction manifest
- [configure](../../configure) - build and HDF5 selection script
- [Doxyfile](../../Doxyfile) - Doxygen configuration
- [.github/run_tests.sh](../../.github/run_tests.sh) - local replay test driver
- [Unit_Tests/sample_table/generate_simple_table.py](../../Unit_Tests/sample_table/generate_simple_table.py) - sample table generator
- [tools/compose/compose_to_grhayl.py](../../tools/compose/compose_to_grhayl.py) - CompOSE EOS converter
- [implementations/GRHayLib/src/make.code.defn](../../implementations/GRHayLib/src/make.code.defn) - GRHayLib thorn manifest

## See Also

- Depends on: [Workflows](../workflows.md#living-source-change)
- See also: [Build And Configure](build-and-configure.md)
- See also: [Public API Map](public-api-map.md)
- See also: [Tests](../tests/index.md)
