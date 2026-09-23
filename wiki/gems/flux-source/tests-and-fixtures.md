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
| Hybrid HLLE flux replay | [Unit_Tests/unit_test_hybrid_flux.c](../../../Unit_Tests/unit_test_hybrid_flux.c) | `hybrid_flux_input.bin`, `hybrid_flux_output.bin`, `hybrid_flux_output_pert.bin` | Direct replay of hybrid and hybrid-entropy variants, all three directions. Also checks wave-bound/error contracts and independent asymmetric fixed-bound HLLE algebra. Speeds are fixture inputs, not recomputed. Stored outputs predate flooring roundoff-negative wave speeds at zero. Samples whose stored `cmin` or `cmax` is negative (roughly half of the per-direction sample comparisons) use a fixed relative/absolute compatibility envelope instead of `ghl_pert_test_fail`: it has no perturbation allowance but a much larger absolute floor (`16*DBL_EPSILON`), so it is usually looser, especially for small fluxes. Nonnegative-speed samples retain perturbation replay. |
| Tabulated HLLE flux replay | [Unit_Tests/unit_test_tabulated_flux.c](../../../Unit_Tests/unit_test_tabulated_flux.c) | `tabulated_flux_input.bin`, `tabulated_flux_output.bin`, `tabulated_flux_output_pert.bin` | Direct replay of tabulated and tabulated-entropy variants, all directions. Also checks wave-bound/error contracts and independent asymmetric fixed-bound HLLE algebra. It initializes LS220, then uses `ghl_test_compute_h_and_cs2`. Stored outputs predate the same wave-speed flooring, and negative-speed samples (roughly half of the per-direction sample comparisons) use the same kind of envelope. Its relative tolerance is much looser than the hybrid one, so tabulated negative-speed samples accept relative differences several orders of magnitude larger than `ghl_pert_test_fail` typically allows. |
| ET Legacy flux/source replay | [Unit_Tests/unit_test_ET_Legacy_flux_source.c](../../../Unit_Tests/unit_test_ET_Legacy_flux_source.c) | `ET_Legacy_flux_source_input.bin`, `ET_Legacy_flux_source_output.bin`, `ET_Legacy_flux_source_output_pert.bin` | Direct combined replay: all speed directions, non-entropy hybrid fluxes, flux divergence, and source terms. It installs one test-local combined enthalpy/sound-speed callback. |
| CompOSE production-EOS integration | [Unit_Tests/unit_test_tabulated_eos_compose.c](../../../Unit_Tests/unit_test_tabulated_eos_compose.c) | Caller-supplied HDF5 path; Ubuntu-GCC builds `compose-test-table.h5` in the workspace root and passes it to the test. | Analytic asymmetric magnetized characteristic-speed expectations in all directions, production tabulated table-bound/mutation, identical-state HLLE/entropy flux, and source-term checks. This focused execution exists only in the Ubuntu-GCC `compose-regularized-eos` job. |

Related contract routes:

- [Characteristic speeds contract](characteristic-speeds-contract.md)
- [HLLE flux variant matrix](hlle-flux-variant-matrix.md)
- [Source-term contract](source-terms-contract.md)
- [Flux_Source hub](../flux-source.md)

## Data Generators

| Generator | Files it writes | Fixture role |
| --- | --- | --- |
| [Unit_Tests/data_gen/unit_test_data_hybrid_flux.c](../../../Unit_Tests/data_gen/unit_test_data_hybrid_flux.c) | `hybrid_flux_input.bin`, `hybrid_flux_output.bin`, `hybrid_flux_output_pert.bin` | Generates hybrid face data, calls all speed and direct HLLE variants, then writes outputs. Producer and consumer use the same kernels, so outputs are regression references, not an independent oracle. |
| [Unit_Tests/data_gen/unit_test_data_tabulated_flux.c](../../../Unit_Tests/data_gen/unit_test_data_tabulated_flux.c) | `tabulated_flux_input.bin`, `tabulated_flux_output.bin`, `tabulated_flux_output_pert.bin` | Initializes LS220, replaces production combined enthalpy/sound-speed dispatch with `ghl_test_compute_h_and_cs2`, calls all speed and direct tabulated HLLE variants, and writes same-kernel regression references. |
| [Unit_Tests/data_gen/unit_test_data_ET_Legacy_flux_source.c](../../../Unit_Tests/data_gen/unit_test_data_ET_Legacy_flux_source.c) | `ET_Legacy_flux_source_input.bin`, `ET_Legacy_flux_source_input_pert.bin` | Generates local ET Legacy input-side fixtures for metric, curvature, primitives, and face states. The test consumes downloaded trusted output and perturbed-output fixtures. |

The hybrid and tabulated generators perturb metric, primitive, and stored
wave-speed inputs; the ET Legacy generator perturbs primitive and face-state
inputs only. Each perturbed input gets an independent uniform relative factor
of up to `1e-14` (hybrid, ET Legacy) or `1e-12` (tabulated).
The repository records no conditioning rationale for the difference. Replay
through `ghl_pert_test_fail` scales each family's bar with its own
trusted-versus-perturbed separation, so tabulated replay is correspondingly
looser. Exact cutoffs and envelopes live in `GRHayL/include/ghl_unit_tests.h`
and the replay sources.

Generator source presence does not mean CI regenerates trusted data. The
ordinary runner downloads trusted fixtures before running tests, and workflow
test jobs also download fixture files before execution.

## Characteristic-Speed Evidence

Characteristic speeds are fixture evidence here. The hybrid and tabulated data
generators call every `ghl_calculate_characteristic_speed_dirn*` kernel and
store `cxmin/cxmax`, `cymin/cymax`, and `czmin/czmax` in the flux input
fixtures. The hybrid and tabulated tests consume those stored speed arrays while
checking HLLE output. The ET Legacy flux/source test instead computes `cmin` and
`cmax` at replay time before calling the hybrid HLLE flux functions.

Neither ordinary flux replay recomputes characteristic speeds. ET Legacy does,
but against its test-local EOS callback. The CompOSE integration test separately
covers production tabulated speed dispatch and table-bound/mutation behavior.
`unit_test_hybrid_flux` injects errors at both callback positions in all three
checked characteristic-speed directions and checks that both outputs stay
unchanged. No committed check injects a failure into the production tabulated
callback during a characteristic-speed call.

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
- Workflow flux jobs download LS220 only for the `tabulated_flux` matrix entry.
- [configure](../../../configure) documents `--disable-hdf5` as disabling
  tabulated EOS. It retains the algebraic tabulated HLLE source files, while
  filtering table-dependent tests and data generators, including the
  `tabulated_flux` test and generator.

## Runner And Workflow Status

- [.github/run_tests.sh](../../../.github/run_tests.sh) runs `./configure -r`,
  `make tests datagen`, sets `LD_LIBRARY_PATH`, downloads fixture files, then runs
  `./test/unit_test_ET_Legacy_flux_source`, `./test/unit_test_hybrid_flux`, and
  `./test/unit_test_tabulated_flux`.
- The runner downloads Flux fixtures from the immutable TestData reference
  recorded in `.github/et-legacy-testdata-ref`. Its exit trap removes only
  run-created paths, including its private expected-error work directory, and
  preserves preexisting paths.
- The workflow files under [.github/workflows/](../../../.github/workflows/)
  have a `flux` job matrix for `hybrid_flux` and `tabulated_flux`. The shared
  compile action builds `make tests datagen`, but test execution still uses
  downloaded fixtures rather than regenerated trusted data.
- ET Legacy workflow matrices include `flux_source`; that is
  `unit_test_ET_Legacy_flux_source`, separate from the workflow matrix entry
  named `HLL_flux`.

Every compiler workflow configures these jobs. Push and pull-request triggers
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

- Magnetized HLLE terms have replay coverage but no independent analytic
  oracle; focused asymmetric HLLE checks use zero magnetic field.
- ET Legacy flux/source replay uses equal spacings (`dx = dy = dz = 0.1`);
  direction-specific inverse spacing is checked only by the selector
  self-check, not by an anisotropic replay.
- Production tabulated mutation is checked for characteristic speeds, not for
  an HLLE call.
- Legacy generic `ghl_calculate_HLLE_fluxes_dirn0/1/2` globals have storage but
  no Core assignment or repository Unit Test call.
- No focused test checks that successful HLLE or source-term calls leave
  non-output conservative fields unchanged; failure-path tests check that all
  fields stay unchanged.
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
