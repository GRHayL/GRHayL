# Neutrino M1 approximation limits and use cases

The current Radiation gem is a deliberately bounded grey transport model. Its
contracts support useful local neutrino number/energy/momentum operations, but
they do not establish a spectrally complete neutrino-transport calculation or
a complete downstream evolution framework. The current capability boundary is
defined by [`ghl_m1.h`](../../../GRHayL/include/ghl_m1.h#L845-L965), the
[provider header](../../../GRHayL/include/ghl_neutrino_rate_provider.h#L14-L27),
and the [integration contract](../../../GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md#current-evidence-boundary).

## Approximation hierarchy

The current model makes these reductions:

1. **Grey, one-group:** each species has one number density, one energy
   density, and one energy-flux vector. There is no group index or retained
   spectrum.
2. **Three species:** `nu_e`, `anti-nu_e`, and an already-summed four-flavor
   `nu_x`. Heavy flavors are not separately evolved.
3. **M1 angular closure:** the E/F moments are closed by an algebraic M1
   pressure tensor. The unresolved angular distribution is not reconstructed.
4. **Frozen-rate local coupling:** the provider supplies a validated rate
   bundle; the pointwise source solve freezes metric, fluid primitives, and
   rates. It does not refresh EOS/rates or call Con2Prim internally.
5. **Grey pair collision:** electron-flavor pair/plasmon/bremsstrahlung fields
   use a first-order occupancy-product number reaction and separate grey energy
   weighting.

These are implementation properties, not optional precision switches. The
host can choose provider channels and source-policy controls, but it cannot
turn the gem into a multigroup or flavor-resolved solver by configuration.

## Defensible current claims

Within the above hierarchy, the implementation can provide:

- local grey transport of number, energy, and momentum for the three modeled
  species;
- provider-defined charged-current, nucleon-scattering, pair,
  bremsstrahlung, and plasmon channel contributions;
- weak-equilibrium-target relaxation when the provider supplies valid
  `n_eq`/`J_eq` and consistent rate identities;
- charged-current electron-lepton and matter `Y_e` bookkeeping; and
- coupled electron-flavor pair number conservation and equal/opposite local
  energy-momentum exchange.

“Provider-defined” matters: the channel mask and backend determine which rates
are active, and the table-free reference backend is a deterministic test model.
The explicit NRPyLeakage initializer selects the table-backed production
provider, subject to its EOS/HDF5 requirements. See
[`ghl_neutrino_rate_provider.h`](../../../GRHayL/include/ghl_neutrino_rate_provider.h#L136-L220)
and [`ghl_neutrino_rate_provider.c`](../../../GRHayL/Radiation/Neutrinos/ghl_neutrino_rate_provider.c#L756-L936).

## Explicit non-claims

The current Radiation M1 surface does not claim to provide:

- multigroup spectral evolution, spectral thermalization, or group-to-group
  energy redistribution;
- separate `nu_mu`, `anti-nu_mu`, `nu_tau`, and `anti-nu_tau` states;
- flavor oscillations or mixing;
- exact inelastic neutrino-electron scattering;
- exact recoil or weak-magnetism corrections;
- complete many-body, heavy-nucleus, or nuclear-composition physics beyond
  what a selected provider actually supplies;
- an exact energy/angle-resolved inverse pair kernel; or
- precision CCSN, PNS, disk, or merger predictions merely because the local
  kernels return finite results.

Some of these topics appear in the neutrino whitepaper as interaction
motivation or future work. They must remain labeled as adaptable context or
roadmap material unless a current provider, source contract, active manifest,
and verification path establish them.

## M1-specific limits

M1 evolves only the retained angular moments. Effects that depend on the full
angular distribution—such as beam crossing, sharp shadows, and detailed
annihilation angular correlations—are therefore outside what can be inferred
from the local M1 state. The whitepaper discussion of these effects is useful
motivation, not a test-backed accuracy guarantee.

The current canonical neutrino transport path is four-point blended Rusanov
with metric light-cone speeds. Older photon-oriented whitepapers describe HLL,
optical-depth speed caps, and a reduced number current; those are not current
neutrino behavior. The current transport distinction is recorded in
[`TRACEABILITY.md`](../../../GRHayL/Radiation/TRACEABILITY.md#claim-boundary),
[`ghl_m1_four_point_blended_rusanov.c`](../../../GRHayL/Radiation/ghl_m1_four_point_blended_rusanov.c#L1-L145),
and [`ghl_m1_neutrino_number_flux.c`](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_number_flux.c#L15-L78).

Optional thick-limit/diffusion helpers exist as separate shared operations;
their presence does not change the canonical four-point neutrino method. They
also do not add spectral or independent number diffusion. See the
[optional diffusion leaf](m1-thick-limit-and-optional-diffusion.md).

## Evidence boundary

The local tests establish selected algebraic and transactional behavior:

- provider channel mapping, multiplicity, equilibrium targets, and rate
  validation in [`unit_test_m1_rate_provider.c`](../../../Unit_Tests/unit_test_m1_rate_provider.c);
- pair number/energy behavior and exchange in
  [`unit_test_m1_neutrino_source_update.c`](../../../Unit_Tests/unit_test_m1_neutrino_source_update.c);
- seeded closure, number-current, source, and matter-coupling invariants in
  [`unit_test_m1_neutrino_seeded_invariants.c`](../../../Unit_Tests/unit_test_m1_neutrino_seeded_invariants.c);
- transport fixtures and current-method comparisons in
  [`unit_test_m1_neutrino_rusanov_flux.c`](../../../Unit_Tests/unit_test_m1_neutrino_rusanov_flux.c)
  and [`unit_test_m1_thcm1_blended_rusanov.c`](../../../Unit_Tests/unit_test_m1_thcm1_blended_rusanov.c).

The [M1 test guide](../../../Unit_Tests/README.m1.md) and
[compatibility evidence leaf](compatibility-evidence.md) state the limits of
that evidence. These tests do not prove a particular host mesh, AMR scheme,
time integrator, downstream thorn, or full-evolution physical validation.

## Historical status

- **Current:** grey three-species state, M1 closure, provider-selected channels,
  frozen local rates, joint grey pair model, and source/exchange boundaries.
- **Adaptable context:** the whitepapers' approximation hierarchy and scenario
  motivation for exploratory grey transport.
- **Superseded:** photon LTE, photon callback microphysics, old HLL-centered
  neutrino transport, and Phase 1's reduced number-current design.
- **Future:** multigroup transport, richer closure/angle models, inelastic and
  recoil physics, weak-magnetism/many-body corrections, and more complete
  downstream verification campaigns.
