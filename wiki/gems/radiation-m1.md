# Radiation M1 Gem

## Purpose

The Radiation M1 gem provides grey, one-group, three-species neutrino transport,
explicit radiation-matter source terms, the canonical four-point blended Rusanov
face operation, local implicit updates, and shared closure, moment, stress,
diagnostic, and diffusion helpers. Generic M1 HLL transport is not part of the
neutrino surface. The public diffusion helpers are optional operations; the
canonical neutrino four-point transport route does not apply a separate
diffusion correction.

The grey three-species neutrino M1 surface is scoped by
[neutrino grey physics scope](radiation-m1/neutrino-grey-physics-scope.md),
governed by the [neutrino reuse strategy](radiation-m1/neutrino-reuse-strategy.md),
and carried forward by the [post-phase-1 roadmap](radiation-m1/post-phase1-roadmap.md).

## Read First

- [Neutrino M1 contract — closure and realizability](radiation-m1/neutrino-m1-contract.md#fixed-numerical-path)
- [Neutrino M1 contract — transport](radiation-m1/neutrino-m1-contract.md#four-point-transport)
- [Neutrino M1 contract — source update](radiation-m1/neutrino-m1-contract.md#transactional-source-update)
- [Neutrino M1 contract — failure and fallback](radiation-m1/neutrino-m1-contract.md#failure-and-fallback-behavior)
- [Neutrino M1 contract — transport and rates](radiation-m1/neutrino-m1-contract.md)
- [Neutrino rate provider](radiation-m1/rate-provider-contract.md)
- [Grey electron-flavor pair collision model](../../GRHayL/Radiation/PAIR_SOURCE_MODEL.md)
- [Neutrino grey physics scope](radiation-m1/neutrino-grey-physics-scope.md)
- [Neutrino reuse strategy](radiation-m1/neutrino-reuse-strategy.md)
- [Post-phase-1 roadmap](radiation-m1/post-phase1-roadmap.md)
- [`GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md`](../../GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md) for
  the transport policy, validation boundary, and configuration/validation
  route
- [`GRHayL/Radiation/TRACEABILITY.md`](../../GRHayL/Radiation/TRACEABILITY.md)
  for the per-family implementation map
- [M1 tests and fixtures](radiation-m1/tests-and-fixtures.md) for the nine
  scoped executables, retained fixture families, provider fixture, CI routes,
  and coverage limits
- [`docs/raw/Radiation.dox`](../../docs/raw/Radiation.dox) for the Doxygen
  group and public behavior summary
- [M1 unit-test guide](../../Unit_Tests/README.m1.md) for the scoped runner,
  frozen THC_M1 comparisons, provenance, and coverage inventory
- [`GRHayL/include/ghl_m1.h`](../../GRHayL/include/ghl_m1.h)
- [`GRHayL/include/ghl_neutrino_rate_provider.h`](../../GRHayL/include/ghl_neutrino_rate_provider.h)
- [Neutrinos hub](neutrinos.md) for the NRPyLeakage raw-rate bridge and
  radiation structs
- [Evolution Equation Map](../physics/evolution-equation-map.md)

## Physics, Mathematics, And Numerical Methods

The contract pages above define the public boundary. These focused leaves
explain the equations and algorithms behind that boundary, with current source
and test references:

### Physics

- [Species, state, and scope](radiation-m1/neutrino-species-state-and-scope.md)
- [Interaction channels](radiation-m1/neutrino-interaction-channels.md)
- [Equilibrium and rate semantics](radiation-m1/neutrino-equilibrium-and-rate-semantics.md)
- [Pair source model](radiation-m1/neutrino-pair-source-model.md)
- [Lepton and matter exchange](radiation-m1/neutrino-lepton-and-matter-exchange.md)
- [Approximation limits and use cases](radiation-m1/neutrino-approximation-limits-and-use-cases.md)

### Mathematics

- [3+1 radiation equations](radiation-m1/m1-3plus1-radiation-equations.md)
- [Comoving moments and source projection](radiation-m1/m1-comoving-moments-and-source-projection.md)
- [Four-dimensional Minerbo closure](radiation-m1/m1-four-dimensional-minerbo-closure.md)
- [Realizability repair](radiation-m1/m1-realizability-repair.md)
- [Number current and transport](radiation-m1/m1-number-current-and-transport.md)
- [Neutrino source equations](radiation-m1/m1-neutrino-source-equations.md)
- [Exchange equations](radiation-m1/m1-exchange-equations.md)

### Numerical methods

- [Finite-volume and face flux](radiation-m1/m1-finite-volume-and-face-flux.md)
- [Four-point blended Rusanov](radiation-m1/m1-four-point-blended-rusanov.md)
- [Geometry and wave-speed policy](radiation-m1/m1-geometry-and-wave-speed-policy.md)
- [Thick limit and optional diffusion](radiation-m1/m1-thick-limit-and-optional-diffusion.md)
- [Implicit Newton and failure policy](radiation-m1/m1-implicit-newton-and-failure-policy.md)
- [Source-update branches and rollback](radiation-m1/m1-source-update-branches-and-rollback.md)
- [Host stages and volume-weighted integration](radiation-m1/m1-host-stage-and-volume-weighted-integration.md)
- [Diagnostics and verification](radiation-m1/m1-diagnostics-and-verification.md)

## Public Headers

- [`GRHayL/include/ghl_m1.h`](../../GRHayL/include/ghl_m1.h)
- [`GRHayL/include/ghl_neutrino_rate_provider.h`](../../GRHayL/include/ghl_neutrino_rate_provider.h)
- [`GRHayL/include/ghl_radiation.h`](../../GRHayL/include/ghl_radiation.h)
- [`GRHayL/include/ghl_nrpyleakage.h`](../../GRHayL/include/ghl_nrpyleakage.h)
- [`GRHayL/include/ghl.h`](../../GRHayL/include/ghl.h)

Both new headers are in the install list in
`GRHayL/include/make.code.defn`.

Key public surface (the installed headers remain the complete declaration
surface):

- Shared M1 kernels used by the neutrino path cover closure, realizability
  repair, comoving moments, metric wave speeds, physical/Rusanov fluxes,
  stress-energy, diagnostics, Newton/projection, `Jthick`, optional diffusion,
  and both pointwise and volume-weighted prepared-face transport support.
- Neutrino transport: `ghl_m1_compute_neutrino_number_flux`,
  `ghl_m1_compute_neutrino_physical_number_flux`,
  `ghl_m1_compute_neutrino_rusanov_flux`,
  `ghl_m1_compute_neutrino_four_point_transport_flux`,
  `ghl_m1_compute_neutrino_four_point_volume_weighted_transport_flux`,
  `ghl_m1_compute_neutrino_diffusion_flux`,
  `ghl_m1_validate_neutrino_rates`, and
  `ghl_m1_repair_neutrino_state`
- Neutrino source/update: `ghl_m1_compute_neutrino_interaction_sources`,
  `ghl_m1_update_neutrino_number_backward_euler`,
  `ghl_m1_try_neutrino_explicit_thin_update`,
  `ghl_m1_solve_neutrino_source_update`,
  `ghl_m1_solve_neutrino_implicit_homogeneous_update`,
  `ghl_m1_solve_neutrino_pair_source_update`,
  `ghl_m1_compute_neutrino_lepton_increment`,
  `ghl_m1_neutrino_compute_implicit_residual`,
  `ghl_m1_neutrino_compute_implicit_jacobian`,
  `ghl_m1_neutrino_build_trial_state`, and
  `ghl_m1_neutrino_check_trial_admissibility`
- Rate provider: `ghl_neutrino_rate_provider_initialize_default`,
  `ghl_neutrino_rate_provider_cache_initialize`,
  `ghl_neutrino_rate_provider_initialize_nrpyleakage`,
  `ghl_neutrino_rate_provider_compute_cell`
- Private NRPyLeakage rate adapter:
  `ghl_m1_nrpyleakage_compute_raw_rates_from_thermo`

The raw-rate bridge is private to the rate provider's Radiation build; it is
not part of `ghl_nrpyleakage.h` and does not extend the legacy scalar leakage
API.
Neutrino species keep the existing leakage naming `nue`, `anue`, `nux`
(`ghl_m1_neutrino_species_count = 3`).

## Implementation Paths

- Shared M1 kernels consumed by neutrino code: `GRHayL/Radiation/`
- Neutrino M1: `GRHayL/Radiation/Neutrinos/`
- Rate provider: `GRHayL/Radiation/Neutrinos/ghl_neutrino_rate_provider.c`
- Build manifest: `GRHayL/Radiation/make.code.defn`,
  `GRHayL/Radiation/Neutrinos/make.code.defn`
- Test runner and fixtures: `Unit_Tests/run_m1_tests.sh`,
  `Unit_Tests/data/m1_thcm1/`, and `.github/actions/run_m1/action.yml`

## Test Paths

- Scoped runner: `Unit_Tests/run_m1_tests.sh`
- Test sources: `Unit_Tests/unit_test_m1_*.c` and
  `Unit_Tests/unit_test_rusanov_flux.c`
- Stored fixture consumers: `Unit_Tests/data/m1_thcm1/`
- Test guide and fixture-package validation: `Unit_Tests/README.m1.md` and
  `Unit_Tests/data/m1_thcm1/audit_package.py`
- Normal runner and dedicated action: `.github/run_tests.sh` and
  `.github/actions/run_m1/action.yml`

## Key Contracts

- Native transport uses the four-point blended Rusanov face operation. It
  applies the metric light-cone speed to `{N, E, Fx, Fy, Fz}` and densitizes
  the final face flux once with the face `sqrt_detgamma`.
- The Induction and Flux_Source HLL/HLLE routines are unrelated to the
  neutrino M1 transport surface.
- The canonical neutrino transport operation has no separate diffusion
  correction or optical-depth wavespeed cap. Optional public `Jthick` and
  diffusion helpers remain available and are tested separately.
- The rate provider computes one frozen rate bundle per grey
  neutrino species. Both built-in backends accept only
  `nu_x_multiplicity == 4` and return an already-summed heavy-lepton bundle.
  The raw-rate bridge exposes one single-species `nux` channel; the provider
  applies the multiplicity once when building the already-summed bundle.
  Hosts must not multiply the returned heavy-lepton state or exchange again.
- Electron pair, plasmon, and bremsstrahlung emissivities remain separate
  from charged-current coefficients and require the paired source update.
  Single-species source operations reject these channels.
- In a `--disable-hdf5` build, tabulated-EOS rate recovery returns
  `ghl_error_used_disabled_hdf5` without publishing output.
- The neutrino implicit Jacobian uses deterministic finite differences of the
  frozen-rate residual. It does not call Con2Prim, recompute opacities, or
  promise an analytic Jacobian.

## Common Edit Routes

- Change a neutrino transport, source, or implicit routine: update the
  relevant shared `GRHayL/Radiation/` or neutrino
  `GRHayL/Radiation/Neutrinos/` source, `GRHayL/include/ghl_m1.h`,
  `GRHayL/Radiation/make.code.defn`,
  the contract, and the traceability map together.
- Add a header or directory: add it to the module and install manifests.
- Change the rate-provider boundary: update
  `ghl_neutrino_rate_provider.h`, `ghl_nrpyleakage.h`, and the provider tests.
- Change an M1 test, fixture schema, or CI selection: update
  [M1 tests and fixtures](radiation-m1/tests-and-fixtures.md),
  `Unit_Tests/README.m1.md`, `Unit_Tests/run_m1_tests.sh`, and the matching
  fixture metadata or CI action.

## Drift Risks

- Generic M1 HLL transport is outside the neutrino API. Separate diffusion is
  a public optional helper, but it is not part of the canonical neutrino
  four-point transport route.
- Unit tests establish unit and invariant behavior only; they are not
  cross-code continuum or Cactus host-integration evidence.

## Do Not Duplicate

Keep equations and per-function documentation in the installed headers and
source. Keep the complete API/validation boundary in the integration contract,
the implementation map in the traceability document, and test/fixture details
in the tests-and-fixtures leaf. This page is for edit routing.
