# Neutrino Grey Physics Scope

## Routing Purpose

Route questions about what the grey, one-group, three-species neutrino M1 model
can and cannot represent to this local scope summary. Use this page when asking
about the interaction inventory, the emissivity/absorption asymmetry, weak
equilibrium and `Y_e` sign rules, the approximation hierarchy, capability verdicts
for BNS/PNS/disk/CCSN, omitted interactions, and which claims are defensible
versus non-claims. This page is the scientific-scope authority for the neutrino
M1 gem; it does not specify file-level API or build truth.

## Public API

The physics whitepaper names no new API. It scopes the physical meaning of the
already-installed surface (declared in [`ghl_m1.h`](../../../GRHayL/include/ghl_m1.h)
and
[`ghl_neutrino_rate_provider.h`](../../../GRHayL/include/ghl_neutrino_rate_provider.h)):

- State: grey three-species `N/E/F_i` per species.
- Frozen rates: `ghl_m1_neutrino_rates` consumed by the source kernels.
- Number source: `eta_N - kappa_a_N * N / Gamma_N`.
- Charged-current `dYe_matter` exchange.
- Provider channels: charged current, nucleon scattering, pair,
  bremsstrahlung, plasmon.

The scope page constrains what those quantities are allowed to mean.

## Caller Contract

- The model is permanently grey, one-group, three-species. A caller must not
  expect group-indexed or spectral behavior from this gem.
- The provider supplies frozen, Kirchhoff-consistent bundles; the Radiation
  kernels do not own production weak-rate physics.
- `Y_e` bookkeeping follows the charged-current subset. Total radiation-number
  exchange remains a separate diagnostic and is not folded into `Y_e`.
- Any caller claim about a physics result must be bounded to what this page
  marks defensible.

## Evidence Status

**Authority label: repository-local scientific scope.** This page records the
approximation boundary used by the current provider and source contracts. Source
and API facts still defer to the installed headers and build lists.

Scope and capability summary (distilled, not copied):

- Representable: grey energy, number, and flux per species; charged-current
  electron-lepton exchange; the interaction channels the provider exposes;
  equilibrium/relaxation behavior at the grey level.
- Not representable: spectral (group) structure, flavor mixing, inelastic
  redistribution, full many-body corrections, heavy-nucleus composition effects
  beyond the provider's model, and the non-equilibrium neutrino field that
  processes like inelastic nu-e scattering require.
- Emissivity/absorption asymmetry and the weak-equilibrium / `Y_e` sign rules are
  the load-bearing approximations a result inherits; a claim must state which
  of these it relies on.
- Capability verdicts are scenario-bounded: the model is scoped for the
  grey three-species regime; BNS/PNS/disk/CCSN statements are defensible only
  within the approximation hierarchy this page records.

## Evidence Links

- [Neutrino M1 contract](neutrino-m1-contract.md) — local neutrino API contract.
- [Rate-provider contract](rate-provider-contract.md) — the channel/physics
  inventory the scope page bounds.
- [Radiation M1 gem](../radiation-m1.md) — hub.

## Ground Truth References

- [`GRHayL/include/ghl_m1.h`](../../../GRHayL/include/ghl_m1.h)
- [`GRHayL/include/ghl_neutrino_rate_provider.h`](../../../GRHayL/include/ghl_neutrino_rate_provider.h)
- [`GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_sources.c`](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_sources.c)
- [`GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md`](../../../GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md)
