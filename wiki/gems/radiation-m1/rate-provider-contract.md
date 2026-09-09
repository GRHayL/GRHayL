# Radiation M1 Rate Provider Contract

Purpose: route questions about the grey neutrino rate-provider boundary. This
page covers the public provider API in
[`ghl_neutrino_rate_provider.h`](../../../GRHayL/include/ghl_neutrino_rate_provider.h)
and the landed staged provider implementation in
[`ghl_neutrino_rate_provider.c`](../../../GRHayL/Radiation/Neutrinos/ghl_neutrino_rate_provider.c).

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

The public provider API exposes:

- `ghl_neutrino_rate_provider_context`: configuration, channel mask,
  failure/table-bound policies, unit scales, per-channel scales,
  `nu_x_multiplicity`, and `use_tabulated_eos`.
- `ghl_neutrino_rate_provider_cache`: thermodynamic state, EOS-derived
  quantities, and cached per-species rates.
- `ghl_neutrino_rate_provider_diagnostics`: failures, table-bound hits,
  cache hits/misses, clamped inputs, recovery counts, active channel mask, and
  last error and last recovery status, plus last-call beta/Kirchhoff relative
  mismatch values and per-species validity flags.
- `ghl_neutrino_rate_provider_initialize_default`.
- `ghl_neutrino_rate_provider_cache_initialize`.
- `ghl_neutrino_rate_provider_initialize_nrpyleakage`.
- `ghl_neutrino_rate_provider_compute_cell`.

The default initializer selects `ghl_neutrino_rate_backend_reference`, keeps
`use_tabulated_eos` false in every build configuration, and is usable without
an EOS table. The explicit
`ghl_neutrino_rate_provider_initialize_nrpyleakage` initializer selects
`ghl_neutrino_rate_backend_nrpyleakage`, requires HDF5 and a Stellarcollapse
tabulated EOS, enables every implemented channel, and defaults both policies to
abort. Both built-in backends require exactly `nu_x_multiplicity = 4` and
return an already-summed `nu_x` bundle.

For electron-flavor transport, scalar `kappa_tr` excludes the separated pair
channels. The host adds their partner-dependent inverse-energy opacity before
constructing the face opacity, as specified in the
[pair collision model](../../../GRHayL/Radiation/PAIR_SOURCE_MODEL.md#opacity-supplied-to-transport).

## Channels

`channel_mask` can enable:

- `ghl_neutrino_rate_channel_charged_current`
- `ghl_neutrino_rate_channel_nucleon_scattering`
- `ghl_neutrino_rate_channel_pair`
- `ghl_neutrino_rate_channel_bremsstrahlung`
- `ghl_neutrino_rate_channel_plasmon`

The reference provider maps enabled channels into its deterministic grey
bundle. The production NRPyLeakage backend maps exposed Ruffert channels:
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
[grey pair collision model](../../../GRHayL/Radiation/PAIR_SOURCE_MODEL.md)
for the equations, time splitting, and energy/number weighting.
Charged-current lepton bookkeeping comes from the provider's charged-current
number fields and the species lepton weights.

The raw beta emissivity is retained as a provider diagnostic, not added
independently to the M1 source. Because the Ruffert beta-emission and
charged-current absorption fits use different blocking averages, their local
equilibrium Kirchhoff reconstructions are checked diagnostically within a
factor of ten (`|eta_beta-eta_cc|/max(eta_beta,eta_cc) <= 0.9`) at the sample
EOS point. The public bundle remains exactly Kirchhoff-consistent through
`eta_N_cc = kappa_a_N_cc*n_eq`; the two approximations are never averaged into
a new formula. This sample-point comparison checks internal rate consistency;
it is not physical validation of the weak-rate model.

On successful table-backed provider calls the runtime diagnostic is
`abs(beta-kirchhoff)/max(abs(beta),abs(kirchhoff),DBL_MIN)`. Only electron
flavors are valid; `nu_x` and reference-backend entries remain invalid. The
unmasked charged-current reconstruction makes the diagnostic independent of
the active thermal or charged-current mask.

Channel toggles are deterministic and are part of provider state. Radiation
does not inspect the channel mask; it only sees the final frozen rates.

## Table Bounds And Table-Free Mode

When `use_tabulated_eos` is true, the provider reads table bounds from
`ghl_eos_parameters`, may recover temperature through tabulated EOS helpers,
and computes thermodynamic quantities through current `NRPyEOS_*` entry
points. If matter inputs fall outside table bounds, the provider either aborts
or clamps according to `table_bounds_policy`. Clamping increments diagnostics.
The reference backend preserves its historical pressure/energy/chemical-
potential lookup followed by the composition lookup. NRPyLeakage alone uses
the single six-quantity chemical-potential/composition interpolation.

When `use_tabulated_eos` is false, only the reference provider uses deterministic
primitive-state formulas for controlled tests and table-free paths. This mode
is the default and is useful for unit tests and verifiers; it is not a
production-calibrated weak-rate table. The production backend must be selected
through the explicit NRPyLeakage initializer.

## Failure Policy

Provider failure behavior is explicit:

- `ghl_neutrino_rate_failure_abort`: return the underlying failure.
- `ghl_neutrino_rate_failure_transparent`: return transparent rates.
- `ghl_neutrino_rate_failure_equilibrium`: return an explicitly degraded,
  Kirchhoff-consistent relaxation bundle. Number and energy absorption use the
  positive caller-configured `equilibrium_recovery_rate` in inverse code-time
  units, emissivities target tiny positive `n_eq`/`J_eq`, electron-flavor
  number coupling is attributed to charged-current exchange, and scattering
  remains zero. This is a safety recovery, not calibrated microphysics.
- `ghl_neutrino_rate_failure_hold_last`: return the most recent bundle only
  when its primitive key, full provider snapshot, EOS pointer, and
  host-managed EOS generation match. Failures before a valid recovered
  `(rho,T,Ye)` key exists are not eligible for hold-last.

Non-abort recoveries increment provider diagnostics and set `last_recovery`;
`last_error` retains the underlying recovered failure. Hold-last recovery depends
on a valid provider cache. Transparent and equilibrium recoveries still return
bundles that must pass Radiation's rate validation before transactional
publication. Invalid fallback arithmetic returns the underlying failure and
leaves the caller's rate array unchanged.

## Cache And Diagnostics

The provider cache is caller-owned and must be initialized with
`ghl_neutrino_rate_provider_cache_initialize` before first use. Its layout intentionally changed to carry
separate thermodynamic and final-rate keys plus a complete provider snapshot
and EOS pointer. The context carries the host-managed EOS generation. Cache
validity and values are committed transactionally only after all species pass
rate validation; exact final-rate hits and `hold_last` require matching
recovered `(rho, T, Ye)`, context, EOS pointer, and generation. Hosts must
advance `eos_generation` after in-place EOS/table mutation.

Provider diagnostics record:

- total failures and last error;
- table-bound hits and clamped inputs;
- cache hits and misses;
- transparent, equilibrium, and hold-last recoveries;
- active channel mask and the most recent call's explicit recovery status;
- beta/Kirchhoff relative mismatches and validity flags from the last fully
  validated table-backed provider result.

These diagnostics are provider-owned. They complement, but do not replace,
`ghl_m1_neutrino_diagnostics.provider_validation_failures` from the Radiation
source/validation path.

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
includes recovered `(rho,T,Ye)`, every context field including backend and
mask, EOS pointer, and `eos_generation`. Failed calls publish neither partial
rates, partial cache state, nor partial beta/Kirchhoff diagnostics. Invalid
backend/table configuration and disabled
HDF5 are boundary errors and are not hidden by a recovery policy.

## `nu_x_multiplicity`

`nu_x` is a lumped heavy-lepton species. The provider context records
`nu_x_multiplicity`, which must be exactly `4` for both built-in backends. The
value denotes the summed heavy-lepton content for `nu_mu`, `anti-nu_mu`,
`nu_tau`, and `anti-nu_tau`; it is not a caller-selectable per-species scale.

The production backend applies the factor exactly once to extensive equilibrium
and emission quantities. It does not multiply opacities or mean energy. The
reference backend likewise returns the already-summed bundle. Radiation kernels
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
- [`GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md`](../../../GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md)
