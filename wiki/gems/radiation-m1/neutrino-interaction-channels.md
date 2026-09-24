# Neutrino interaction channels

The provider exposes five channel selectors, but the channel selector is not
itself a source equation. It controls which provider physics contributes to a
frozen rate bundle; the Radiation kernels then consume that bundle. The public
channel enum and ownership rule are in
[`ghl_neutrino_rate_provider.h`](../../../GRHayL/include/ghl_neutrino_rate_provider.h).

## Channel boundary

The current channel set is:

| channel | current rate representation | source/composition consequence |
| --- | --- | --- |
| charged current | electron-flavor number and energy absorption/emission; a separate number subset `kappa_a_N_cc`, `eta_N_cc` | changes electron lepton number and therefore can change `Y_e` |
| nucleon scattering | scalar `kappa_s`, described as isoenergetic scattering opacity | contributes to transport extinction; it is not a number-reaction or charged-current source |
| pair | process-indexed electron-flavor `eta_N_pair` and `eta_E_pair` | requires the joint `nu_e`/`anti-nu_e` source operation |
| plasmon | process-indexed pair emissivity fields for electron flavors; aggregate scalar contribution for `nu_x` | follows the same paired electron-flavor path when represented in process arrays |
| bremsstrahlung | process-indexed pair emissivity fields for electron flavors; aggregate scalar contribution for `nu_x` | follows the same paired electron-flavor path when represented in process arrays |

The table describes the current provider/source contract, not a guarantee that
every backend call enables every channel. `channel_mask` is provider state;
Radiation sees only the final validated rates. The production mapping is in
[`ghl_neutrino_rate_provider.c`](../../../GRHayL/Radiation/Neutrinos/ghl_neutrino_rate_provider.c),
and the deterministic reference mapping is in the same file
([`compute_staged_rates`](../../../GRHayL/Radiation/Neutrinos/ghl_neutrino_rate_provider.c)).

## Charged current

Charged-current rates are available for `nu_e` and `anti-nu_e`. The provider
stores the number subset separately:

```text
eta_N_cc = kappa_a_N_cc * n_eq.
```

This subset is used only for electron-lepton bookkeeping. It must not be
reconstructed from the total radiation-number increment when pair channels are
active. The energy absorption coefficient contributes to the ordinary grey
energy source; the number absorption coefficient contributes to the number
source. The current source equations are implemented in
[`ghl_m1_neutrino_sources.c`](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_sources.c)
and the subset is validated in
[`ghl_m1_neutrino_rates.c`](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_rates.c).

The implementation does not expose a separate detailed reaction state for
each beta process. The provider owns the production weak-rate formulas and
maps them into the grey bundle. A raw beta-emissivity consistency diagnostic
in the table-backed provider is not an additional M1 source term.

## Nucleon scattering

The provider exposes `kappa_s` as isoenergetic scattering opacity and forms the
transport coefficient

```text
kappa_tr = kappa_a_E + kappa_s.
```

Scattering damps the comoving flux in the grey collision operator without
creating or destroying radiation number. In a moving Eulerian frame its
projection can appear in the energy and momentum source components, while the
comoving scattering energy exchange remains zero. The current projection is

```text
Q       = eta_E - kappa_a_E J
S_E     = Q W + kappa_tr H_n
S_i     = Q W V_i - kappa_tr H_i.
```

Here `J`, `H_i`, and `H_n` are the M1 comoving moments; this is a source
projection, not a claim that scattering changes electron fraction. The source
implementation is authoritative
([`ghl_m1_neutrino_sources.c`](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_sources.c)).

## Pair, plasmon, and bremsstrahlung emission

For electron flavors, these channels are retained as process-indexed emission
data rather than being silently folded into the scalar one-species source:

```text
eta_N_pair[c], eta_E_pair[c],
c in {pair, plasmon, bremsstrahlung}.
```

The number emissivity for each process must match between `nu_e` and
`anti-nu_e`; the energy emissivities may differ. Their inverse reaction depends
on the partner occupancy, so the host must include the partner-dependent
inverse-energy opacity when preparing an electron-flavor transport face. The
scalar `kappa_tr` in the validated electron bundle is not overwritten with
that partner-dependent quantity. These rules are specified in
[`PAIR_SOURCE_MODEL.md`](../../../GRHayL/Radiation/PAIR_SOURCE_MODEL.md#provider-coefficients)
and have separate ownership. The pair kernel enforces the shared
number-emissivity equality while performing the coupled pair update. The
provider/rate-validation boundary owns the validity of the supplied rate
bundle, and the host owns preparation of the partner-dependent transport
opacity and preservation of the validated bundle. The pair kernel does not
inspect or enforce those caller-prepared face quantities; its equality check
is in
[`ghl_m1_neutrino_pair_source.c`](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_pair_source.c).

For the lumped heavy flavor, the provider keeps the process arrays zero and
maps its already-summed pair/plasmon/bremsstrahlung content into aggregate
scalar emissivity/absorption coefficients. The heavy-flavor multiplicity is
applied exactly once. The current public field semantics are documented in
[`ghl_m1.h`](../../../GRHayL/include/ghl_m1.h), and the production raw
mapping is visible in
[`ghl_m1_nrpyleakage_kernel.c`](../../../GRHayL/Radiation/Neutrinos/ghl_m1_nrpyleakage_kernel.c).

These channels are grey emission/absorption representations. The provider
does not deliver an energy-resolved pair spectrum to the M1 source.

## Validation identities

For the scalar one-body bundle, the current validator requires finite,
nonnegative fields and the binary64-consistent identities

```text
eta_N    = kappa_a_N    * n_eq
eta_E    = kappa_a_E    * J_eq
eta_N_cc = kappa_a_N_cc * n_eq
kappa_tr = kappa_a_E + kappa_s
J_eq     = n_eq * mean_energy.
```

Representational underflow is recorded diagnostically when the correctly
rounded product is zero; it is not permission to invent a nonzero rate. The
validator also rejects pair fields in a one-species electron source call and
rejects nonzero pair fields for `nu_x`. See
[`ghl_m1_neutrino_rates.c`](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_rates.c).

## Historical whitepaper status

- **Current:** the five channel names, provider ownership, separate electron
  pair fields, aggregate `nu_x`, and source-level distinction between
  charged-current number exchange and scattering.
- **Adaptable context:** the neutrino interaction whitepaper's inventory of
  free-nucleon charged current, nucleon scattering, pair annihilation,
  bremsstrahlung, and plasmon processes.
- **Superseded:** any Phase 1 design that treated all electron-flavor number
  emission as a one-species aggregate source when separated pair fields are
  present.
- **Future/non-claim:** inelastic redistribution, exact spectral pair
  kinetics, weak-magnetism/recoil corrections, and many-body or heavy-nucleus
  enhancements are not implied by this channel enum.

## Evidence

The channel mapping and identities are exercised by
[`unit_test_m1_rate_provider.c`](../../../Unit_Tests/unit_test_m1_rate_provider.c).
The coupled pair and composition behavior is exercised by
[`unit_test_m1_neutrino_source_update.c`](../../../Unit_Tests/unit_test_m1_neutrino_source_update.c).
