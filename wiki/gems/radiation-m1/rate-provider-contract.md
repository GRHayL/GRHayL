# Radiation M1 Rate Provider Contract

Purpose: route questions about the grey neutrino rate-provider boundary. This
page covers the public provider API in
[`ghl_neutrino_rate_provider.h`](../../../GRHayL/include/ghl_neutrino_rate_provider.h)
and the provider orchestration in
[`ghl_neutrino_rate_provider.c`](../../../GRHayL/Radiation/Neutrinos/ghl_neutrino_rate_provider.c).
The private [rate backend](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_rate_backend.c)
owns the tabulated-EOS callbacks and physical rate assembly, with file-scope
enabled and disabled-HDF5 implementations.

## Boundary

Radiation consumes frozen `ghl_m1_neutrino_rates` bundles. The provider owns
the rate source, EOS/table access, unit conversion, channel policy, table-bound
policy, failure recovery, cache state, and provider diagnostics.

The pointwise Radiation kernels do not own production weak-rate physics. They
validate and consume frozen bundles through the neutrino M1 kernels, returning
`ghl_error_m1_microphysics_failure` when a bundle violates the contract.

The design rationale for keeping the frozen-rate solve separate from the
pointwise kernels — the reuse rubric, file-by-file disposition, and the
provider design contract — is distilled in
[Neutrino reuse strategy](neutrino-reuse-strategy.md). The scientific scope the
provider's channel set is bounded by lives in
[Neutrino grey physics scope](neutrino-grey-physics-scope.md).

## Public Provider Surface

The installed API exposes a production context, caller-owned cache and
diagnostics, the default and explicit NRPyLeakage initializers, cache
initialization, and `ghl_neutrino_rate_provider_compute_cell`.
The context owns the channel mask, failure/table-bound policies,
`nu_x_multiplicity`, EOS generation and equilibrium recovery opacity.
The production provider always uses a tabulated EOS and returns
`ghl_error_used_disabled_hdf5` without publishing rates in a no-HDF5 build.
Synthetic rates and reference-only scales belong to test-local support,
not the installed API or libghl.

Both initializers select the same production provider. The default failure
and table-bound policies return the original error. `nu_x_multiplicity` is
exactly four and the returned heavy-flavor bundle is already summed.

For electron-flavor transport, scalar `kappa_tr` excludes the separated pair
channels. The host adds their partner-dependent inverse-energy opacity before
constructing the face opacity, as specified in the
[pair collision model](../../../docs/raw/Radiation_pair_source_model.md#opacity-supplied-to-transport).

## Channels

`channel_mask` can enable:

- `ghl_neutrino_rate_channel_charged_current`
- `ghl_neutrino_rate_channel_nucleon_scattering`
- `ghl_neutrino_rate_channel_pair`
- `ghl_neutrino_rate_channel_bremsstrahlung`
- `ghl_neutrino_rate_channel_plasmon`

The production provider maps exposed Ruffert channels:
charged-current channels contribute to electron flavors, nucleon scattering
contributes to all species. Electron-flavor pair, plasmon, and bremsstrahlung
emissivities are retained separately in `eta_N_pair` and `eta_E_pair`; scalar
electron-flavor absorption/emission contains only independent charged-current
channels. The lumped heavy species retains aggregate detailed-balance
absorption. `kappa_tr` is formed in the exact order `kappa_a_E + kappa_s`.
Number-weighted scattering is diagnostic raw data only and never enters
number absorption.

The independent number source is `eta_N-kappa_a_N*N/Gamma_N`. Electron-flavor
pair reactions require the coupled operation and its shared reaction extent;
one-species source calls reject those rates. See the
[grey pair collision model](../../../docs/raw/Radiation_pair_source_model.md)
for the equations, time splitting, and energy/number weighting.
Charged-current lepton bookkeeping comes from the provider's charged-current
number fields and the species lepton weights.

The raw beta emissivity is retained as a provider diagnostic, not added
independently to the M1 source. On successful table-backed provider calls,
the provider stores
`abs(beta-kirchhoff)/max(abs(beta),abs(kirchhoff),DBL_MIN)` in
`beta_kirchhoff_relative_mismatch` and sets
`beta_kirchhoff_mismatch_valid` for electron flavors. The `nu_x` entry remains invalid. The public bundle remains exactly
Kirchhoff-consistent through `eta_N_cc = kappa_a_N_cc*n_eq`; the beta and
Kirchhoff approximations are never averaged into a new formula. The
unmasked charged-current reconstruction makes the diagnostic independent of
the active thermal or charged-current mask.

Channel toggles are deterministic and are part of provider state. Radiation
does not inspect the channel mask; it only sees the final frozen rates.
The production raw kernel evaluates only enabled channel rates, while it keeps
the equilibrium Fermi moments and electron-flavor beta/Kirchhoff consistency
diagnostics independent of the mask. Consequently, an empty mask publishes
finite equilibrium targets with zero interaction coefficients and pair-process
emissivities. A disabled channel's unused Fermi integral cannot fail the
provider call.

## Table bounds

The provider reads bounds from `ghl_eos_parameters`, recovers temperature when
required, and obtains composition and chemical potentials through the
initialized tabulated-EOS callbacks. The bounds policy either returns the
underlying error or clamps inputs, recording the clamp in diagnostics.

## Failure and physical recovery

`ghl_neutrino_rate_failure_return_error` preserves caller rates on failure.
The transparent and equilibrium policies can publish a recovery only after
successful same-cell thermodynamic reconstruction and positive, finite
physical Fermi equilibrium moments. The transparent bundle has zero
interaction coefficients. The equilibrium bundle uses the configured
`equilibrium_recovery_rate` and those physical number/energy targets, with
Kirchhoff-consistent emissivities and charged-current electron-flavor number
exchange. No fabricated mean energy or tiny equilibrium target is substituted.

Every recovery still returns the original error and records `last_recovery`.
A host accepting a recovery must check that diagnostic. If thermodynamic or
physical-target reconstruction fails, rates remain unchanged. Disabled HDF5
and invalid configuration are boundary errors. The former hold-last policy
was removed because an exact validated cache hit already returns success
before any failing rate calculation.

## Cache and diagnostics

Initialize each cache with `ghl_neutrino_rate_provider_cache_initialize`.
Exact reuse includes recovered `(rho,T,Ye)`, the complete context, EOS pointer
and host-managed `eos_generation`. Advance that generation after in-place EOS
mutation. Cache values commit only after every species validates.

Diagnostics record failures, bounds/clamp events, cache hits/misses,
transparent/equilibrium recoveries, active channels, last error/recovery and
beta/Kirchhoff mismatch values with validity flags. These complement
Radiation's bundle-validation diagnostics. Cache and diagnostics are mutable
caller-owned records; concurrent callers require separate records or external
synchronization.

## Production units and equilibrium policy

The NRPyLeakage backend uses literal neutrino number. Number density converts
with `L0^3`; MeV energy density with `C_MeV L0^3/E0`; mean energy with
`C_MeV/E0`; opacity with `L0`; and the corresponding emissivities include one
factor of `t0`. Here `L0`, `t0`, and `E0=M0 c^2` come from
`ghl_nrpyleakage.h`. Its equilibrium degeneracies are local and tau-free:
`(mu_e-muhat)/T`, its negative, and zero for `nu_e`, `anti-nu_e`, and `nu_x`.
Leakage diffusion suppression, source terms, and luminosities never enter an M1
bundle.

## Cache ownership and failure

Cache and diagnostics records are caller-owned and mutable. A host must provide
one record per cell or thread, or externally synchronize access. Exact reuse
includes recovered `(rho,T,Ye)`, every context field including the
mask, EOS pointer, and `eos_generation`. Failed calls publish neither partial
rates, partial cache state, nor partial beta/Kirchhoff diagnostics. Invalid
provider/table configuration and disabled
HDF5 are boundary errors and are not hidden by a recovery policy.

## `nu_x_multiplicity`

`nu_x` is a lumped heavy-lepton species. The provider context records
`nu_x_multiplicity`, which must be exactly `4`. The
value denotes the summed heavy-lepton content for `nu_mu`, `anti-nu_mu`,
`nu_tau`, and `anti-nu_tau`; it is not a caller-selectable per-species scale.

The production backend applies the factor exactly once to extensive equilibrium
and emission quantities. It does not multiply opacities or mean energy. Radiation kernels
and hosts do not consume the multiplicity field or reapply the factor.

## Physics inventory and limitations

NRPyLeakage supplies its existing Ruffert charged-current, free-nucleon
scattering, pair, bremsstrahlung, and plasmon algebra. Omitted physics includes
heavy-nucleus capture/scattering, weak magnetism, recoil, inelastic
neutrino-electron scattering, many-body corrections, and spectral
thermalization. Electron-flavor pair number stoichiometry is enforced by the
paired source operation, using the documented grey occupancy-product
approximation. Spectral/angular pair kernels and a coupled matter/radiation
implicit solve remain outside this collision model.

## Local Ground Truth

- [`GRHayL/include/ghl_neutrino_rate_provider.h`](../../../GRHayL/include/ghl_neutrino_rate_provider.h)
- [`GRHayL/include/ghl_m1.h`](../../../GRHayL/include/ghl_m1.h)
- [`GRHayL/Radiation/Neutrinos/ghl_neutrino_rate_provider.c`](../../../GRHayL/Radiation/Neutrinos/ghl_neutrino_rate_provider.c)
- [`GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_rates.c`](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_rates.c)
- [`GRHayL/Radiation/Neutrinos/make.code.defn`](../../../GRHayL/Radiation/Neutrinos/make.code.defn)
- [`docs/raw/Radiation_integration_contract.md`](../../../docs/raw/Radiation_integration_contract.md)
