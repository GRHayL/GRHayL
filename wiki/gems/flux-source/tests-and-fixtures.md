# Flux_Source Tests And Fixtures

## Routing Purpose

Use this page to route Flux_Source test and fixture questions for hydrodynamic
HLLE fluxes, ET Legacy flux/source replay, characteristic-speed fixture
evidence, and HDF5 table needs. Equations, public API contracts, and generated
formula boundaries stay in the Flux_Source contract pages and source evidence.
This page is a test router, not a replacement for test files.

## Direct Flux_Source Test Coverage

Status vocabulary on this page:

- **direct replay:** test source calls the named routine and compares fixtures;
- **fixture generation:** data-generator source calls it, but normal jobs do not
  regenerate fixtures;
- **runner-listed / CI-configured:** command exists in runner/workflow source,
  not proof of a current execution or pass;
- **compiled-unrun:** requires a reviewer-observed binary; tracked source alone
  cannot establish this status;
- **coverage gap:** no focused evidence for the stated behavior.

| Test | Main evidence | Fixture files read | Notes |
| --- | --- | --- | --- |
| Hybrid HLLE flux replay | [Unit_Tests/unit_test_hybrid_flux.c](../../../Unit_Tests/unit_test_hybrid_flux.c) | `hybrid_flux_input.bin`, `hybrid_flux_output.bin`, `hybrid_flux_output_pert.bin` | Direct replay of hybrid and hybrid-entropy variants, all three directions. Speeds are fixture inputs, not recomputed. |
| Tabulated HLLE flux replay | [Unit_Tests/unit_test_tabulated_flux.c](../../../Unit_Tests/unit_test_tabulated_flux.c) | `tabulated_flux_input.bin`, `tabulated_flux_output.bin`, `tabulated_flux_output_pert.bin` | Direct replay of tabulated and tabulated-entropy variants, all directions. It initializes LS220 and computes face thermodynamics, then replaces production `ghl_compute_h_and_cs2` with `ghl_test_compute_h_and_cs2`. |
| ET Legacy flux/source replay | [Unit_Tests/unit_test_ET_Legacy_flux_source.c](../../../Unit_Tests/unit_test_ET_Legacy_flux_source.c) | `ET_Legacy_flux_source_input.bin`, `ET_Legacy_flux_source_output.bin`, `ET_Legacy_flux_source_output_pert.bin` | Direct combined replay: all speed directions, non-entropy hybrid fluxes, flux divergence, and source terms. It also installs a test-local enthalpy/sound-speed callback. |

Related contract routes:

- [Characteristic speeds contract](characteristic-speeds-contract.md)
- [HLLE flux variant matrix](hlle-flux-variant-matrix.md)
- [Source-term contract](source-terms-contract.md)
- [Flux_Source hub](../flux-source.md)

## Data Generators

| Generator | Files it writes | Fixture role |
| --- | --- | --- |
| [Unit_Tests/data_gen/unit_test_data_hybrid_flux.c](../../../Unit_Tests/data_gen/unit_test_data_hybrid_flux.c) | `hybrid_flux_input.bin`, `hybrid_flux_output.bin`, `hybrid_flux_output_pert.bin` | Generates hybrid face data, calls all speed and direct HLLE variants, then writes outputs. Producer and consumer use the same kernels, so outputs are regression references, not an independent oracle. |
| [Unit_Tests/data_gen/unit_test_data_tabulated_flux.c](../../../Unit_Tests/data_gen/unit_test_data_tabulated_flux.c) | `tabulated_flux_input.bin`, `tabulated_flux_input_pert.bin`, `tabulated_flux_output.bin`, `tabulated_flux_output_pert.bin` | Initializes LS220, replaces production enthalpy/sound-speed dispatch with the test helper, calls all speed and direct tabulated HLLE variants, and writes same-kernel regression references. Replay does not read `tabulated_flux_input_pert.bin`. |
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

Neither ordinary flux replay recomputes characteristic speeds. ET Legacy does,
but against its test-local EOS callback. This leaves production tabulated
dispatch, mutation, and error behavior outside direct speed coverage.

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
  script, then removes root-level `*.bin`, `*.h5`, and `*.bz2` at the end. Run
  it only in a clean disposable worktree, not beside caller-owned root data
  with those suffixes.
- The workflow files under [.github/workflows/](../../../.github/workflows/)
  have a `flux` job matrix for `hybrid_flux` and `tabulated_flux`. The shared
  compile action builds `make tests datagen`, but test execution still uses
  downloaded fixtures rather than regenerated trusted data.
- ET Legacy workflow matrices include `flux_source`; that is
  `unit_test_ET_Legacy_flux_source`, separate from the workflow matrix entry
  named `HLL_flux`.

All five workflow files configure these jobs. Push and pull-request triggers
ignore `wiki/**` and Markdown-only changes; scheduled triggers remain. Thus a
KB-only change does not itself exercise Flux_Source on push/PR.

Treat entries above as source-configured execution paths, not observed
current-worktree passes. Record compile and execution results separately during
verification.

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

## Coverage Gaps

- Generic `ghl_calculate_HLLE_fluxes_dirn0/1/2` globals have storage but no
  assignment, call, or test.
- No focused test checks ignored EOS callback errors, tabulated primitive
  mutation, zero `cmin + cmax`, or conservative fields left untouched.
- No direct test connects Flux_Source characteristic-speed output to Induction
  HLL input.

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
