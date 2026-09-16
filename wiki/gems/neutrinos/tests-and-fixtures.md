# Neutrinos Tests And Fixtures

Purpose: route NRPyLeakage test setup, fixture data, and CI commands back to
repo-local ground truth. Source, tests, `configure`, and CI scripts remain the
authority.

## Shared Harness

Every NRPyLeakage unit-test file includes `Unit_Tests/nrpyleakage_main.h`.
That helper owns the common CLI:

```sh
./test/unit_test_<nrpyleakage test> <EOS table path> <test key>
```

The accepted keys are:

- `0`: call `generate_test_data(&eos)` and write local binary reference data.
- `1`: call `run_unit_test(&eos)` and attempt comparison against existing
  binary fixtures; see the replay assertion strength below.

Any other key is rejected by the shared harness. The harness zero-initializes
the EOS and lets `ghl_initialize_tabulated_eos_functions_and_params(...)`
select the tabulated StellarCollapse implementation with sound-speed cleaning
disabled (`GRHayL/GRHayL_Core/initialize_eos.c`). It passes atmosphere values
for `rho_b`, `Y_e`, and `T`, and frees table memory with
`ghl_tabulated_free_memory(&eos)` after either mode completes.

Ground truth:

- `Unit_Tests/nrpyleakage_main.h`
- `GRHayL/include/ghl_nrpyleakage.h`
- `GRHayL/include/ghl_radiation.h`

## Fixture Pairs

Fixture names are hard-coded by each test and are downloaded by
`.github/run_tests.sh` for normal CI replay:

- `nrpyleakage_optically_thin_gas_unperturbed.bin`
- `nrpyleakage_optically_thin_gas_perturbed.bin`
- `nrpyleakage_constant_density_sphere_unperturbed.bin`
- `nrpyleakage_constant_density_sphere_perturbed.bin`
- `nrpyleakage_luminosities_unperturbed.bin`
- `nrpyleakage_luminosities_perturbed.bin`

The same set can be written locally by running the corresponding executable with
key `0`; CI and standard unit-test replay use key `1`. These producers live in
the test translation units themselves; there are no matching
`Unit_Tests/data_gen/unit_test_data_nrpyleakage_*.c` targets.

The historical published `GRHayL/TestData` copies predate the current
diffusion-time, bremsstrahlung-density, heavy-species, and nucleon-blocking
corrections. The owner accepted regenerated outputs from the corrected model as
the replacement golden baseline. Those replacements are installed in the local
`TestData` checkout and await publication to the separate repository. Until
publication, remote CI still downloads the historical data. No new scenario or
tolerance change accompanies the replacement.

## Blocking-Correction Qualification

Current blocking derives free-nucleon kinetic degeneracies from `rho`, `T`,
`X_n`, and `X_p`. The former `Ye=0.5` composition branch and full-`muhat`
transition quotient no longer exist. Density inversion removes the arbitrary
chemical-energy-zero dependence; stable overlap identities remove the old
full-`muhat` zero-denominator pole and the `Ye=0.5` branch discontinuity.
The shifted charged-current parent kernels are paired to the EOS equilibrium
chemical potential. Those mathematical properties do not establish
interacting-EOS accuracy for the algebraic grey reduction.

`Unit_Tests/unit_test_nrpyleakage_physics.c` supplies persistent, table-free
physics checks. It compares blocking and ordinary beta moments with independent
high-precision references, checks population bounds and the transition
normalization identity, checks both reaction-threshold orientations, and
verifies spectral detailed balance for electron-neutrino and
electron-antineutrino kernels. It also exercises the zero-emission,
finite-normalized Boltzmann, mixed numerator-only underflow, and
paired-subnormal limits of strongly blocked channels, preventing `0/0`,
premature underflow, and overflow. Its
floating-point tolerance follows the published Fermi-fit error bound plus an
allowance for the small fixed operation chain. It runs in no-HDF5 builds and
once per compiler and OS neutrino CI matrix.

During development, the existing luminosity, optically-thin, and
constant-density-sphere key-`0` generation scenarios completed with the new
implementation. Key `0` writes data and is not a fixture replay. The sphere
generator derives a temperature margin from the table lower bound and its
existing perturbation amplitude so its exterior sample stays inside the EOS
domain. The owner accepted the qualification results and authorized regenerated
replay outputs as the golden baseline for the new blocking model.

Independent physical qualification used the beta kernels in
[BNS_NURATES](https://github.com/RelNucAs/bns_nurates) and the six DD2 merger
states published by Chiesa et al., *Phys. Rev. D* **111**, 063053 (2025),
[doi:10.1103/PhysRevD.111.063053](https://doi.org/10.1103/PhysRevD.111.063053).
At the existing dilute SLy4 state, electron-neutrino and
electron-antineutrino number and energy emission differed by `2.0--4.0%`.
After `0.5 s`, `Ye` differed by `0.00559` absolute, specific energy by
`0.0166%`, and temperature by `0.00455%`. RK4 step refinement changed final
`Ye` by `3.4e-12`.

The full DD2 comparison exposes a material dense-EOS limitation. At the
densest point, the common-bare-mass model reconstructed `q=117.10 MeV`, versus
`20.22 MeV` from DD2 mean fields and effective masses. The resulting complete
beta rate and opacity differences were factors of roughly `58--356` for
electron neutrinos and `11--26` for electron antineutrinos. They decrease
toward the lower-density points, but several antineutrino channels still
differ by `12--39%` at the three lowest-density published states. Widened
scratch-only BNS_NURATES quadrature changed the reported full-DD2 channels by
at most `6.4e-6` relative from 96 to 128 nodes.

This record does not invent a universal physical pass limit. Neither project
policy, ILEAS, nor BNS_NURATES supplies a per-kernel or one-zone evolution
tolerance applicable here. ILEAS's scheme-level transport agreement is not such
a bound. After reviewing the measured discrepancies and uncertainties, the
owner accepted them for this approximate leakage model and authorized them as
the golden baseline. This outcome-specific decision retains the dense DD2
limitation. EOS-consistent effective masses or mean-field shifts remain an
optional future accuracy improvement.

All established NRPyLeakage replay fixtures were regenerated from the accepted
model. A second independent generation produced byte-identical files, and the
optically thin, constant-density sphere, and luminosity replays passed against
the installed local `TestData` results.

### Cost Evidence

An isolated alternating benchmark ran 500 process invocations of the existing
optically-thin key-`0` executable for the candidate and an `LD_PRELOAD`
no-blocking stub. Mean wall time was `3.315 s` for the candidate and `3.205 s`
for the stub, a `3.43%` local leakage cost. Using the user-supplied estimate
that leakage consumes `35%` of a BNS run gives the back-of-the-envelope total
runtime impact

$$
0.35\times3.43\%\simeq1.20\%.
$$

The owner accepts this estimate as sufficient performance evidence for the
blocking correction. It supports the requested `1--2%` whole-run budget under
that runtime attribution. It is not an end-to-end BNS measurement. Initialization,
call frequency, compiler, hardware, and the share of leakage work spent in the
three EOS-dependent routines can change the result. The evaluator was selected
for this cost constraint: it uses rational approximations and algebraic
overlaps, with no quadrature, iterative root solve, or extra EOS/table lookup.

## HDF5 And EOS Setup

NRPyLeakage tests require an HDF5-backed tabulated EOS path. The repo-local CI
route downloads `SLy4_3335_rho391_temp163_ye66.h5.bz2`, unpacks it to
`SLy4_3335_rho391_temp163_ye66.h5`, downloads the Neutrinos fixture set,
then runs:

```sh
./test/unit_test_nrpyleakage_physics
./test/unit_test_nrpyleakage_optically_thin_gas SLy4_3335_rho391_temp163_ye66.h5 1
./test/unit_test_nrpyleakage_constant_density_sphere SLy4_3335_rho391_temp163_ye66.h5 1
./test/unit_test_nrpyleakage_luminosities SLy4_3335_rho391_temp163_ye66.h5 1
```

Use `.github/run_tests.sh` as the primary command route for fixture and table
downloads. Do not claim current remote fixture availability from this page; the
repo-local script is the cited route.

No-HDF5 builds set `GHL_DISABLE_HDF5` and apply `configure`'s exact
implementation-source filter. That filter removes many tabulated EOS/Flux_Source
sources but retains selected Con2Prim Tabulated/NN helpers and all NRPyLeakage
implementation files. It also excludes the EOS-table NRPyLeakage tests:

- `Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c`
- `Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c`
- `Unit_Tests/unit_test_nrpyleakage_luminosities.c`

The table-free `Unit_Tests/unit_test_nrpyleakage_physics.c` remains compiled
and runnable without HDF5.

NRPyLeakage implementation sources remain compiled so their early
`ghl_error_used_disabled_hdf5` paths exist, but exclusion of those tests
means those paths have no direct no-HDF5 test. `unit_test_code_error` keys `2`
and `3` cover only invalid Fermi keys and remain available without HDF5.

For wider build context, see `wiki/build-and-ci.md` and `wiki/test-map.md`.

## Test Behavior Map

`Unit_Tests/unit_test_nrpyleakage_physics.c` covers density-derived nucleon
blocking and shifted charged-current algebra without an EOS table or binary
fixture. Independent reference values include an equal-population state, a
near-equal state that exercises cancellation handling, a degenerate state,
an asymmetric dilute-proton state, and all four ordinary beta-moment helper
paths. Exact checks cover population bounds, transition normalization,
zero-shift recovery, threshold orientation, and spectral Kirchhoff pairing for
both beta channels. Cold trace-population states cover common moment underflow,
numerator-only underflow with representable Fermi normalization, the finite
normalized Boltzmann limit, and finite ratios of subnormal moments.

`Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c` covers optically thin
gas source replay. Its RHS calls
`NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms`, divides
`R_source` and `Q_source` by `rho`, and evolves `Y_e` and `eps` through RK4.
Each RK4 substep recomputes temperature through `ghl_tabulated_compute_T_from_eps`.
Generation writes the time series; key `1` replays and calls comparison helpers
against
`nrpyleakage_optically_thin_gas_{unperturbed,perturbed}.bin`.

`Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c` covers opacity and
optical-depth behavior. It computes interior/exterior opacities with
`NRPyLeakage_compute_neutrino_opacities`, fills grid arrays, then iterates
optical-depth updates through `NRPyLeakage_optical_depths_PathOfLeastResistance`
until convergence or iteration limit. Key `1` reads and calls comparison
helpers on
`nrpyleakage_constant_density_sphere_{unperturbed,perturbed}.bin`.

`Unit_Tests/unit_test_nrpyleakage_luminosities.c` has two coverage layers. First,
it checks selected `NRPyLeakage_Fermi_Dirac_integrals` branches for large and
small `z`. Then it replays random metric, primitive, optical-depth, and
luminosity cases through `NRPyLeakage_compute_neutrino_luminosities`, then calls
comparison helpers for `nue`, `anue`, and `nux` values against
`nrpyleakage_luminosities_{unperturbed,perturbed}.bin`.

`Unit_Tests/unit_test_code_error.c` belongs here only for direct invalid
Fermi-Dirac key coverage: it calls `NRPyLeakage_Fermi_Dirac_integrals(-1, ...)`
for both `z < 1e-3` and `z > 1e-3` cases and maps those keys to
`ghl_error_invalid_fermi_dirac_integral_key`.

### Replay Assertion Strength

Every NRPyLeakage replay test consumes the boolean returned by
`ghl_pert_test_fail` and calls `ghl_error` on the first mismatch, reporting the
field, its location, and the trusted, computed, and perturbed values:

- Optically thin replay fails on time, `Y_e`, `eps`, or `T` mismatch, reporting
  the evolution time.
- Constant-density-sphere replay fails on any of the six opacity and six
  optical-depth arrays, reporting the grid coordinates and flat index. Its
  final validation loop is serial so that the first failure is reported.
- Luminosity replay fails on `nue`, `anue`, or `nux` mismatch, reporting the
  row index.

All three pass the unperturbed fixture value as `trusted`, the recomputed value
as `computed`, and the perturbed fixture value as `perturbed`, matching the
helper contract. The optically thin test also checks every
`ghl_tabulated_compute_T_from_eps` and `ghl_tabulated_compute_eps_from_T`
status through `ghl_abort_if_error`, so a failed inversion cannot silently
evaluate a right-hand side at a stale temperature.

The luminosity test's explicit Fermi-Dirac assertion rejects a nonfinite
computed integral or a nonfinite reference in addition to an out-of-tolerance
difference.

Remaining bounded gaps:

- Luminosity explicitly checks only high-branch Fermi keys `0`, `1`, `2` and
  low-branch keys `0`, `1`; valid keys `3` through `5` are indirect only.
- The constant-density-sphere stencils remain flat, so unequal-direction face
  metrics are still untested even though its neighbor pairs now match the
  `NRPyLeakage_optical_depths_PathOfLeastResistance` minus-then-plus signature.
- Combined source-term opacity outputs are computed but never compared in the
  optically thin test.
- No Neutrinos test deliberately exercises EOS-error propagation,
  disabled-HDF5 returns, non-finite final luminosity/source outputs, or the
  `EnsureFinite`/`robust_isfinite` fallback paths.

## CI Routes

`.github/run_tests.sh` is the compact end-to-end route for Neutrinos fixture
downloads and test commands. The concrete workflow example is
`.github/workflows/github-actions-Ubuntu-gcc.yml`: its `neutrinos` job matrix
lists `nrpyleakage_optically_thin_gas`,
`nrpyleakage_constant_density_sphere`, and `nrpyleakage_luminosities`, then
downloads matching `_unperturbed.bin` and `_perturbed.bin` files plus the SLy4
table before running each test with key `1`.

Matching Neutrinos jobs exist in other workflow files. Before expanding line
citations beyond the Ubuntu GCC example, verify exact locations with:

```sh
rg -n "neutrinos|nrpyleakage" .github/workflows
```

The exact workflow set is Ubuntu GCC, Ubuntu Clang, Ubuntu Intel, macOS GCC,
and macOS Clang. Each uses its supported OS versions and every test name.
Coverage upload is active for GCC and Ubuntu Clang Neutrinos jobs; it is
commented out for Intel and macOS Clang. Workflow selection makes these
numerical comparisons run; it does not by itself qualify the downloaded
fixtures.

The aggregate runner deletes working-directory `*.bin`, `*.h5`, and `*.bz2`
files at the end. Individual workflow jobs rely on fresh workspaces and do not
include an explicit Neutrinos cleanup step.

## Repo-Local References

- `Unit_Tests/nrpyleakage_main.h`
- `Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c`
- `Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c`
- `Unit_Tests/unit_test_nrpyleakage_luminosities.c`
- `Unit_Tests/unit_test_code_error.c`
- `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c`
- `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities.c`
- `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_luminosities.c`
- `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_Fermi_Dirac_integrals.c`
- `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_optical_depths_PathOfLeastResistance.c`
- `GRHayL/include/ghl_nrpyleakage.h`
- `GRHayL/include/ghl_radiation.h`
- `configure`
- `.github/run_tests.sh`
- `.github/workflows/github-actions-Ubuntu-gcc.yml`
- `wiki/build-and-ci.md`
- `wiki/test-map.md`
