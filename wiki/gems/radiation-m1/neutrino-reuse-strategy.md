# Neutrino Reuse Strategy

## Routing Purpose

Route questions about the Radiation Infrastructure Reuse Strategy to this local
design summary. Use this page when asking why the neutrino M1 work extends the shared
Radiation framework instead of building neutrino-only code, what the reuse rubric
is, how the file-by-file disposition works, why the frozen-rate solve is a
separate component, and what the provider design contract is. This page is a
design/doctrine reference; the landed boundary is the
[rate-provider contract](rate-provider-contract.md).

## Public API

The reuse strategy names no new API. It explains the ownership split that the
landed surface already embodies:

- Shared Radiation framework: the `ghl_m1_*` pointwise closure, repair,
  moments, geometry, wave-speed, prepared-face, stress-energy, diagnostics,
  and Newton families used by neutrino transport.
- Neutrino-specific kernels: `ghl_m1_neutrino_*` (state, rates, number flux,
  repair, prepared-face transport, sources, implicit, exchange, and lepton
  increment).
- Provider boundary: `ghl_neutrino_rate_provider_*` — the frozen-rate solve
  kept separate from the pointwise Radiation kernels.

## Caller Contract

- All later work extends the shared Radiation framework where the neutrino
  surface can use it; a new neutrino-specific component requires a separate
  contract and justification.
- Reuse existing neutrino-compatible infrastructure directly before adding a
  new neutrino-specific implementation.
- The frozen-rate solve is a separate component so that production weak-rate
  physics never enters the pointwise kernels; Radiation consumes frozen bundles.
- The number-current closure rationale and the repair budget are part of the
  design contract a host and provider must respect.

## Evidence Status

**Authority label: repository-local design/rationale authority.** This is
doctrine, not current API or build authority.
Cross-check the landed ownership boundary against
`ghl_neutrino_rate_provider.h` and the Neutrinos build list before treating a
reused name as current.

Recorded discrepancy (doctrine vs landed):

- Parked and reference-only components are not API authority. The reuse rubric
  is the stable doctrine; the specific file-by-file disposition must be checked
  against the current build lists.
- The landed provider is the concrete realization of the "frozen-rate solve is
  separate" rule; see the
  [rate-provider contract](rate-provider-contract.md) for the public surface and
  the [post-phase-1 roadmap](post-phase1-roadmap.md) for how the doctrine
  governs later phases.

## Evidence Links

- [Rate-provider contract](rate-provider-contract.md) — landed provider boundary.
- [Post-phase-1 roadmap](post-phase1-roadmap.md) — doctrine applied to later phases.
- [Radiation M1 gem](../radiation-m1.md) — hub.

## Ground Truth References

- [`GRHayL/include/ghl_neutrino_rate_provider.h`](../../../GRHayL/include/ghl_neutrino_rate_provider.h)
- [`GRHayL/Radiation/Neutrinos/make.code.defn`](../../../GRHayL/Radiation/Neutrinos/make.code.defn)
- [`GRHayL/Radiation/Neutrinos/ghl_neutrino_rate_provider.c`](../../../GRHayL/Radiation/Neutrinos/ghl_neutrino_rate_provider.c)
- [`GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md`](../../../GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md)
