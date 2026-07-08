# Flux_Source Tests And Fixtures

## Routing Purpose

Use this page to route Flux_Source test and fixture questions for hydrodynamic
HLLE fluxes, ET Legacy flux/source replay, characteristic-speed fixture
evidence, and HDF5 table needs. Equations, public API contracts, and generated
formula boundaries stay in the Flux_Source contract pages and source evidence.
This page is a test router, not a replacement for test files.

## Direct Flux_Source Test Coverage

| Test | Main evidence | Fixture files read | Notes |
| --- | --- | --- | --- |
| Hybrid HLLE flux replay | [Unit_Tests/unit_test_hybrid_flux.c](../../../Unit_Tests/unit_test_hybrid_flux.c) | `hybrid_flux_input.bin`, `hybrid_flux_output.bin`, `hybrid_flux_output_pert.bin` | Covers hybrid and hybrid-entropy `ghl_calculate_HLLE_fluxes_dirn*` variants for all three directions. The input fixture carries metric, right/left primitives, and generated `cmin`/`cmax` arrays. |
| Tabulated HLLE flux replay | [Unit_Tests/unit_test_tabulated_flux.c](../../../Unit_Tests/unit_test_tabulated_flux.c) | `tabulated_flux_input.bin`, `tabulated_flux_output.bin`, `tabulated_flux_output_pert.bin` | Covers tabulated and tabulated-entropy `ghl_calculate_HLLE_fluxes_dirn*` variants for all three directions. It initializes `LS220_234r_136t_50y_analmu_20091212_SVNr26.h5` and computes pressure, eps, and entropy from temperature before replaying fluxes. |
| ET Legacy flux/source replay | [Unit_Tests/unit_test_ET_Legacy_flux_source.c](../../../Unit_Tests/unit_test_ET_Legacy_flux_source.c) | `ET_Legacy_flux_source_input.bin`, `ET_Legacy_flux_source_output.bin`, `ET_Legacy_flux_source_output_pert.bin` | Replays hybrid Flux_Source flux/source behavior against ET Legacy data. It computes characteristic speeds directly, calls hybrid HLLE fluxes, builds flux-divergence RHS terms, then adds `ghl_calculate_source_terms` output. |

Related contract routes:

- [Characteristic speeds contract](characteristic-speeds-contract.md)
- [HLLE flux variant matrix](hlle-flux-variant-matrix.md)
- [Source-term contract](source-terms-contract.md)
- [Flux_Source hub](../flux-source.md)

## Data Generators

| Generator | Files it writes | Fixture role |
| --- | --- | --- |
| [Unit_Tests/data_gen/unit_test_data_hybrid_flux.c](../../../Unit_Tests/data_gen/unit_test_data_hybrid_flux.c) | `hybrid_flux_input.bin`, `hybrid_flux_output.bin`, `hybrid_flux_output_pert.bin` | Generates hybrid face data, calls `ghl_calculate_characteristic_speed_dirn0/1/2`, stores `cmin`/`cmax` in the input fixture, then writes trusted and perturbed HLLE outputs. |
| [Unit_Tests/data_gen/unit_test_data_tabulated_flux.c](../../../Unit_Tests/data_gen/unit_test_data_tabulated_flux.c) | `tabulated_flux_input.bin`, `tabulated_flux_input_pert.bin`, `tabulated_flux_output.bin`, `tabulated_flux_output_pert.bin` | Initializes the LS220 table, generates tabulated face data, computes characteristic speeds, writes the unperturbed input consumed by the test, and writes trusted and perturbed HLLE outputs. |
| [Unit_Tests/data_gen/unit_test_data_ET_Legacy_flux_source.c](../../../Unit_Tests/data_gen/unit_test_data_ET_Legacy_flux_source.c) | `ET_Legacy_flux_source_input.bin`, `ET_Legacy_flux_source_input_pert.bin` | Generates local ET Legacy input-side fixtures for metric, curvature, primitives, and face states. The test consumes downloaded trusted output and perturbed-output fixtures. |

Generator source presence does not mean CI regenerates trusted data. The
ordinary runner downloads trusted fixtures before running tests, and workflow
test jobs also download fixture files before execution.

## Characteristic-Speed Evidence

Characteristic speeds are fixture evidence here. The hybrid and tabulated data
generators call the three `ghl_calculate_characteristic_speed_dirn*` kernels and
store `cxmin/cxmax`, `cymin/cymax`, and `czmin/czmax` in the flux input
fixtures. The hybrid and tabulated tests consume those stored speed arrays while
checking HLLE output. The ET Legacy flux/source test instead computes `cmin` and
`cmax` at replay time before calling the hybrid HLLE flux functions.

Do not duplicate speed equations in this page. Public behavior and formula
routes belong to [characteristic-speeds-contract.md](characteristic-speeds-contract.md),
[docs/raw/Flux_Source.dox](../../../docs/raw/Flux_Source.dox), and
[docs/raw/derivation.md](../../../docs/raw/derivation.md).

## HDF5 And Table Needs

- `unit_test_hybrid_flux` uses hybrid EOS setup and does not require an HDF5
  table.
- `unit_test_tabulated_flux` requires HDF5-enabled tabulated EOS support and
  `LS220_234r_136t_50y_analmu_20091212_SVNr26.h5` in the working directory.
- `.github/run_tests.sh` downloads
  `LS220_234r_136t_50y_analmu_20091212_SVNr26.h5.bz2` from the repo-visible
  stellarcollapse URL and uncompresses it before `unit_test_tabulated_flux`.
- The workflow flux jobs download the same LS220 table for both
  `hybrid_flux` and `tabulated_flux` matrix entries.
- [configure](../../../configure) documents `--disable-hdf5` as disabling
  tabulated EOS and filters out tabulated sources, tests, and data generators,
  including `tabulated_flux`.

## Runner And Workflow Status

- [.github/run_tests.sh](../../../.github/run_tests.sh) runs `./configure -r`,
  `make tests`, sets `LD_LIBRARY_PATH`, downloads fixture files, then runs
  `./test/unit_test_ET_Legacy_flux_source`, `./test/unit_test_hybrid_flux`, and
  `./test/unit_test_tabulated_flux`.
- The same runner downloads `flux_source/hybrid_flux_*` and
  `flux_source/tabulated_flux_*` from the GRHayL TestData path recorded in the
  script, then removes root-level `*.bin`, `*.h5`, and `*.bz2` at the end.
- The workflow files under [.github/workflows/](../../../.github/workflows/)
  have a `flux` job matrix for `hybrid_flux` and `tabulated_flux`. The shared
  compile action builds `make tests datagen`, but test execution still uses
  downloaded fixtures rather than regenerated trusted data.
- ET Legacy workflow matrices include `flux_source`; that is
  `unit_test_ET_Legacy_flux_source`, separate from the workflow matrix entry
  named `HLL_flux`.

## Not Flux_Source HLLE Coverage

`Unit_Tests/unit_test_HLL_flux.c` and
`Unit_Tests/unit_test_ET_Legacy_HLL_flux.c` are Induction vector-potential HLL
evidence, not Flux_Source HLLE coverage. They route through
[Induction HLL flux contract](../induction/hll-flux-contract.md), where the
public evidence is `ghl_HLL_flux_with_B`, `ghl_HLL_flux_with_Btilde`, and
`ghl_HLL_vars`.

Any future Flux_Source hub or test-map route that treats
`unit_test_HLL_flux.c` as Flux_Source coverage should be corrected to
Induction vector-potential HLL routing. Flux_Source HLLE coverage is the
hydrodynamic `unit_test_hybrid_flux`, `unit_test_tabulated_flux`, and
`unit_test_ET_Legacy_flux_source` set above.

## Evidence Links

- [wiki/gems/flux-source.md](../flux-source.md)
- [wiki/gems/flux-source/characteristic-speeds-contract.md](characteristic-speeds-contract.md)
- [wiki/gems/flux-source/hlle-flux-variant-matrix.md](hlle-flux-variant-matrix.md)
- [wiki/gems/flux-source/source-terms-contract.md](source-terms-contract.md)
- [wiki/gems/induction/hll-flux-contract.md](../induction/hll-flux-contract.md)
- [Unit_Tests/unit_test_hybrid_flux.c](../../../Unit_Tests/unit_test_hybrid_flux.c)
- [Unit_Tests/unit_test_tabulated_flux.c](../../../Unit_Tests/unit_test_tabulated_flux.c)
- [Unit_Tests/unit_test_ET_Legacy_flux_source.c](../../../Unit_Tests/unit_test_ET_Legacy_flux_source.c)
- [Unit_Tests/data_gen/unit_test_data_hybrid_flux.c](../../../Unit_Tests/data_gen/unit_test_data_hybrid_flux.c)
- [Unit_Tests/data_gen/unit_test_data_tabulated_flux.c](../../../Unit_Tests/data_gen/unit_test_data_tabulated_flux.c)
- [Unit_Tests/data_gen/unit_test_data_ET_Legacy_flux_source.c](../../../Unit_Tests/data_gen/unit_test_data_ET_Legacy_flux_source.c)
- [.github/run_tests.sh](../../../.github/run_tests.sh)
- [.github/workflows/](../../../.github/workflows/)
- [configure](../../../configure)
