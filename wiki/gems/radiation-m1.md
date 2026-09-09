# Radiation M1 Gem

## Purpose

The Radiation M1 gem provides grey, one-group, three-species neutrino transport,
explicit radiation-matter source terms, the canonical four-point blended Rusanov
face operation, and local implicit updates. Generic M1 HLL transport and a
separate diffusion correction are not part of the neutrino surface.

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
- [`GRHayL/include/ghl_m1.h`](../../GRHayL/include/ghl_m1.h)
- [`GRHayL/include/ghl_neutrino_rate_provider.h`](../../GRHayL/include/ghl_neutrino_rate_provider.h)
- [Neutrinos hub](neutrinos.md) for the NRPyLeakage raw-rate bridge and
  radiation structs
- [Evolution Equation Map](../physics/evolution-equation-map.md)

## Public Headers

- [`GRHayL/include/ghl_m1.h`](../../GRHayL/include/ghl_m1.h)
- [`GRHayL/include/ghl_neutrino_rate_provider.h`](../../GRHayL/include/ghl_neutrino_rate_provider.h)
- [`GRHayL/include/ghl_radiation.h`](../../GRHayL/include/ghl_radiation.h)
- [`GRHayL/include/ghl_nrpyleakage.h`](../../GRHayL/include/ghl_nrpyleakage.h)
- [`GRHayL/include/ghl.h`](../../GRHayL/include/ghl.h)

Both new headers are in the install list in
`GRHayL/include/make.code.defn`.

Key public surface:

- Shared M1 kernels used by the neutrino path cover closure, realizability
  repair, comoving moments, metric wave speeds, stress-energy, diagnostics,
  and the prepared-face transport support.
- Neutrino: `ghl_m1_compute_neutrino_number_flux`,
  `ghl_m1_compute_neutrino_rusanov_flux`,
  `ghl_m1_compute_neutrino_interaction_sources`,
  `ghl_m1_update_neutrino_number_backward_euler`,
  `ghl_m1_solve_neutrino_implicit_homogeneous_update`,
  `ghl_m1_solve_neutrino_pair_source_update`,
  `ghl_m1_neutrino_compute_implicit_residual`,
  `ghl_m1_neutrino_compute_implicit_jacobian`,
  `ghl_m1_compute_neutrino_lepton_increment`
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

## Key Contracts

- Native transport uses the four-point blended Rusanov face operation. It
  applies the metric light-cone speed to `{N, E, Fx, Fy, Fz}` and densitizes
  the final face flux once with the face `sqrt_detgamma`.
- The Induction and Flux_Source HLL/HLLE routines are unrelated to the
  neutrino M1 transport surface.
- The canonical neutrino transport operation has no separate diffusion
  correction or optical-depth wavespeed cap.
- The rate provider computes one frozen rate bundle per grey
  neutrino species. Both built-in backends accept only
  `nu_x_multiplicity == 4` and return an already-summed heavy-lepton bundle.
  The raw-rate bridge output is one single-species `nux` channel; callers that
  evolve four heavy-lepton species own the aggregation multiplicity.
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

## Drift Risks

- Generic M1 HLL transport and separate diffusion correction are outside the
  neutrino API; do not describe them as supported neutrino methods.
- Unit tests establish unit and invariant behavior only; they are not
  cross-code continuum or Cactus host-integration evidence.

## Do Not Duplicate

Keep equations and per-function documentation in the installed headers and
source. Keep the complete API/validation boundary in the integration contract
and the implementation map in the traceability document. This page is for
edit routing.
