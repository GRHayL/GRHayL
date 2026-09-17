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

The published `GRHayL/TestData` copies contain regenerated outputs for the
current diffusion-time, bremsstrahlung-density, heavy-species, and
nucleon-blocking corrections. CI downloads those replacement goldens. Replay
scenarios remain the established three, but sphere and luminosity comparisons
now use `1024*DBL_EPSILON` with zero absolute tolerance, and the sphere exterior
temperature respects the EOS-table lower bound. These implementation-generated
baselines detect drift; they do not independently validate corrected physics.

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
physics checks. It compares blocking with independent high-precision
references and ordinary beta moments with high-precision evaluations of the
same fitted formulas, checks population bounds and the transition
normalization identity, checks both reaction-threshold orientations, and
verifies spectral detailed balance for electron-neutrino and
electron-antineutrino kernels. Direct analytic identities catch loss of the
`6/c` diffusion factor, the bremsstrahlung `rho^2` scaling, and omission of the
four heavy-lepton species from the matter-energy source. An independent
public-API thick-limit reconstruction uses public opacity, equilibrium energy
density, and a zero-depth free rate to catch wrong species or units at both
production emission call sites. Exact binary64 patterns
exercise finite values, infinities, quiet and signaling NaNs, both signed
zeros, and both signed minimum subnormals. The representation-based checks
avoid dependence on optimizer finite-value assumptions.
It also exercises the zero-emission,
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

These measurements are external qualification evidence. Their scratch drivers,
DD2 state data, and modified BNS_NURATES quadrature are not present in this
checkout, so this repository cannot reproduce those numerical comparisons.

This record does not invent a universal physical pass limit. Neither project
policy, ILEAS, nor BNS_NURATES supplies a per-kernel or one-zone evolution
tolerance applicable here. ILEAS's scheme-level transport agreement is not such
a bound. After reviewing the measured discrepancies and uncertainties, the
owner accepted them for this approximate leakage model and authorized them as
the golden baseline. This outcome-specific decision retains the dense DD2
limitation. EOS-consistent effective masses or mean-field shifts remain an
optional future accuracy improvement.

All established NRPyLeakage replay fixtures were regenerated from the accepted
model. A second generation using the same implementation produced byte-identical files, and the
optically thin, constant-density sphere, and luminosity replays passed against
the installed local `TestData` results.

### Cost Evidence

An earlier alternating benchmark ran 500 whole-process invocations of the
optically-thin key-`0` executable for the candidate and an `LD_PRELOAD`
no-blocking stub. Mean wall times were `3.315 s` and `3.205 s`. Each invocation
also loaded an approximately 879 MB EOS table, so the observed `3.43%`
difference is not a leakage-kernel measurement. Multiplying it by the
user-supplied `35%` leakage share does not establish whole-BNS overhead.

No repeated production-kernel benchmark or end-to-end BNS timing is recorded
in this checkout. Whole-BNS impact therefore remains unmeasured. The evaluator
was selected to control cost: it uses rational approximations and algebraic
overlaps, with no quadrature, iterative root solve, or extra EOS/table lookup.
That design fact supports a low-cost expectation but is not a quantitative
runtime acceptance result.

## HDF5 And EOS Setup

NRPyLeakage tests require an HDF5-backed tabulated EOS path. The repo-local CI
route downloads `SLy4_3335_rho391_temp163_ye66.h5.bz2`, unpacks it to
`SLy4_3335_rho391_temp163_ye66.h5`, downloads the Neutrinos fixture set,
then runs:

```sh
./test/unit_test_nrpyleakage_physics
./test/unit_test_nrpyleakage_classifier_fallback
./test/unit_test_nrpyleakage_optically_thin_gas SLy4_3335_rho391_temp163_ye66.h5 1
./test/unit_test_nrpyleakage_constant_density_sphere SLy4_3335_rho391_temp163_ye66.h5 1
./test/unit_test_nrpyleakage_luminosities SLy4_3335_rho391_temp163_ye66.h5 1
```

Use `.github/run_tests.sh` as the local aggregate command route for fixture and
table downloads. Current workflows do not invoke that script; their neutrino
jobs run the same table-free and table-backed executables directly. Command
presence proves configured intent, while actual execution requires a local log
or workflow result. The script's named TestData paths are the reproducible
publication route for the current fixtures described above.

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
`ghl_error_used_disabled_hdf5` paths exist. The table-free physics test calls
all three public leakage routines in its no-HDF5 build and checks those return
codes directly. `unit_test_code_error` keys `2` and `3` separately cover the
invalid Fermi keys.

For wider build context, see `wiki/build-and-ci.md` and `wiki/test-map.md`.

## Test Behavior Map

`Unit_Tests/unit_test_nrpyleakage_physics.c` covers density-derived nucleon
blocking and shifted charged-current algebra without an EOS table or binary
fixture. Independent blocking reference values include an equal-population
state, a near-equal state that exercises cancellation handling, a degenerate
state, and an asymmetric dilute-proton state. High-precision arithmetic
references cover all four ordinary beta-moment helper paths. Exact checks cover
the stable order-zero Fermi expression, roundoff-sized fraction normalization,
population bounds, transition normalization, zero-shift recovery, threshold
orientation, and spectral Kirchhoff pairing for both beta channels. A
hand-calculated asymmetric stencil checks optical-depth neighbor ordering with
unequal face metrics. Cold trace-population states cover common moment underflow,
numerator-only underflow with representable Fermi normalization, the finite
normalized Boltzmann limit, and finite ratios of subnormal moments.

`Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c` covers optically thin
gas source replay. Its RHS calls
`NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms`, divides
`R_source` and `Q_source` by `rho`, and evolves `Y_e` and `eps` through RK4.
Each RK4 substep recomputes temperature through `ghl_tabulated_compute_T_from_eps`.
Failure injection verifies the four lookup sites return before 0, 1, 2, or 3
right-hand-side calls. A separate injection verifies the initial energy lookup
occurs before opening the fixture output, so failure cannot truncate it.
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
Generation retains each base value, applies a bounded perturbation, and asserts
the pair distance and EOS/input bounds before evaluating the perturbed row.

`Unit_Tests/unit_test_code_error.c` belongs here only for direct invalid
Fermi-Dirac key coverage: it calls `NRPyLeakage_Fermi_Dirac_integrals(-1, ...)`
for both `z < 1e-3` and `z > 1e-3` cases and maps those keys to
`ghl_error_invalid_fermi_dirac_integral_key`.

### Replay Assertion Strength

Every NRPyLeakage replay test calls `ghl_error` on the first comparison mismatch,
reporting the field, its location, and the trusted, computed, and perturbed
values. The optically thin replay uses `ghl_pert_test_fail`. The luminosity and
constant-density-sphere replays use `ghl_pert_test_fail_with_tolerance` with a
local relative floor of `1024*DBL_EPSILON`. Both local policies set the
absolute tolerance to zero, so very small opacity, optical-depth, and
luminosity values remain subject to a relative comparison instead of passing
through the generic `1e-30` absolute floor:

- When both trusted and perturbed fixtures are exactly zero, the local policy
  requires the computed value to be exactly zero. A focused negative control
  verifies that even `DBL_MIN` is rejected in this case.

- Optically thin replay fails on time, `Y_e`, `eps`, or `T` mismatch, reporting
  the evolution time.
- Constant-density-sphere replay fails on any of the six opacity and six
  optical-depth arrays, reporting the grid coordinates and flat index. Its
  final validation loop is serial so that the first failure is reported.
- Luminosity replay fails on `nue`, `anue`, or `nux` mismatch, reporting the
  row index.

All three pass the unperturbed fixture value as `trusted`, the recomputed value
as `computed`, and the perturbed fixture value as `perturbed`, matching the
helper contract. The optically thin test propagates every
`ghl_tabulated_compute_T_from_eps` and `ghl_tabulated_compute_eps_from_T`
status, so a failed inversion cannot evaluate a right-hand side at a stale
temperature or truncate a fixture before the initial lookup succeeds.

The luminosity test's explicit Fermi-Dirac assertion rejects a nonfinite
computed integral or a nonfinite reference in addition to an out-of-tolerance
difference.

The table-free physics test checks valid Fermi keys `0` through `5` on both
sides of the fit branch, a hand-computed stencil with distinct neighbors and
unequal face metrics, and equality of standalone and combined opacity outputs.
It also injects EOS success and failure callbacks, checks the opacity floor and
neutral rate/source fallback for every public output, verifies public
finite-depth heavy-lepton suppression against manufactured regression values,
reconstructs the public heavy-lepton diffusion limit independently, checks
both single-species endpoint orientations against the tiny-positive analytic
limit through all three public APIs, checks the effective-rate
endpoints and bremsstrahlung energy conversion, and
verifies all three disabled-HDF5 returns.

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
