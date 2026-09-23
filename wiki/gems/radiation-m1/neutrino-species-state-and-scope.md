# Neutrino species, state, and scope

This leaf explains what the current Radiation M1 implementation evolves. It is
the physics companion to the [neutrino M1 contract](neutrino-m1-contract.md)
and the [rate-provider contract](rate-provider-contract.md); it does not
replace either API contract.

## Current model

The public model is grey, one-group, and three-species. The enumeration and
state definition are in [`ghl_m1.h`](../../../GRHayL/include/ghl_m1.h#L845-L895):

| code species | physical content | electron-lepton weight |
| --- | --- | ---: |
| `ghl_m1_neutrino_nue` | electron neutrino, `nu_e` | `+1` |
| `ghl_m1_neutrino_anue` | electron antineutrino, `anti-nu_e` | `-1` |
| `ghl_m1_neutrino_nux` | one lumped heavy-lepton species | `0` |

`nu_x` represents the sum of `nu_mu`, `anti-nu_mu`, `nu_tau`, and
`anti-nu_tau`. Both built-in provider backends require
`nu_x_multiplicity == 4` and return that state already summed. A host must not
multiply it again. This multiplicity is a species-model convention, not a
fourth evolved state or a caller-selectable flavor resolution; see the
[provider header](../../../GRHayL/include/ghl_neutrino_rate_provider.h#L14-L27)
and its [provider implementation](../../../GRHayL/Radiation/Neutrinos/ghl_neutrino_rate_provider.c#L528-L549).

For each species, the undensitized state is

```text
U_s = { N_s, E_s, F_{x,s}, F_{y,s}, F_{z,s} }.
```

Here `N` is Eulerian radiation number density, `E` is Eulerian radiation
energy density, and `F_i` is the covariant Eulerian energy flux. The host may
store densitized variables, but pointwise Radiation calls receive undensitized
states; the header documents the required `sqrt(det(gamma))` conversion
([`ghl_m1.h`](../../../GRHayL/include/ghl_m1.h#L882-L895)). The fixed transport
component order is also recorded by the
[neutrino Rusanov implementation](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_rusanov_flux.c#L1-L80).

`N` is an independently transported grey moment. It is not reconstructed by
dividing `E` by a fixed mean energy. Its transport current is derived from the
M1 E/F closure and comoving moments, with an explicit `N == 0` branch, in
[`ghl_m1_neutrino_number_flux.c`](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_number_flux.c#L15-L78).

## What the state can represent

The state and provider bundle support, at grey level:

- energy and momentum transport through `E`, `F_i`, and the M1 pressure
  closure;
- a separate radiation number current through `N`;
- per-species equilibrium targets and absorption/emission relaxation;
- charged-current electron-lepton exchange for the two electron flavors; and
- a joint electron-flavor pair reaction with equal number increments and
  separate energy weighting.

The source kernels consume provider-supplied rates. They do not derive weak
rates from the EOS themselves. The provider owns EOS/table lookup, unit
conversion, weak-equilibrium targets, and channel-specific microphysics; the
boundary is specified in the
[provider header](../../../GRHayL/include/ghl_neutrino_rate_provider.h#L14-L27)
and [integration contract](../../../GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md#ownership-rules).

The library is pointwise and host-neutral. The downstream host owns the mesh,
cell storage, reconstruction, flux divergence, Runge--Kutta or other stage
schedule, boundary conditions, matter recovery, and final publication. The
current source boundary is summarized in the
[M1 integration contract](../../../GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md#ownership-rules).

## What is not a hidden part of the model

There is no group-indexed transport API. A grey `N` and a grey `E` do not retain
an energy spectrum, and the M1 closure does not retain the full angular
distribution. The three-species state also does not resolve heavy-lepton
flavors separately or evolve flavor mixing.

The current production surface is determined by the active shared and
neutrino manifests, not by every file present in the directory. See the
[shared manifest](../../../GRHayL/Radiation/make.code.defn),
[neutrino manifest](../../../GRHayL/Radiation/Neutrinos/make.code.defn), and
[traceability map](../../../GRHayL/Radiation/TRACEABILITY.md#production-boundary).

## Historical whitepaper status

- **Current:** three grey species, `N/E/F_i`, provider-owned frozen rates, and
  the distinction between electron-lepton and heavy-lepton content.
- **Adaptable context:** the neutrino interaction whitepaper's motivation for
  evolving number as well as energy, and its warning that grey M1 is an
  approximation. Those papers remain explanatory context, not API authority.
- **Superseded:** Phase 1's reduced number-current construction and its
  proposed APIs. The current number current reuses the M1 closure and
  comoving moments; it is not simply `N F^i/F`.
- **Future/non-claim:** multigroup spectra, four separately evolved heavy
  flavors, flavor oscillations, and spectral angular transport.

The [source-level traceability map](../../../GRHayL/Radiation/TRACEABILITY.md)
is the final check when a whitepaper statement and this summary disagree.

## Evidence

The focused tests exercise the current state/rate boundary and selected local
physics operations:

- [`unit_test_m1_rate_provider.c`](../../../Unit_Tests/unit_test_m1_rate_provider.c)
  checks species indexing, multiplicity, channel mapping, and validated rate
  bundles.
- [`unit_test_m1_neutrino_source_update.c`](../../../Unit_Tests/unit_test_m1_neutrino_source_update.c)
  checks local source, pair, exchange, and transactional behavior.
- [`unit_test_m1_neutrino_seeded_invariants.c`](../../../Unit_Tests/unit_test_m1_neutrino_seeded_invariants.c)
  checks generated local states and source/current invariants.

These tests are evidence for the operations they call; they do not establish
full mesh evolution or downstream framework coverage.
