# Inelastic Neutrino-Electron/Positron Scattering (Future Work)

Status: future "to do" work. This page is not a record of implemented or
validated functionality. No inelastic neutrino-electron channel is present in
the current checked-in NRPyLeakage source, and none is claimed here.

## Why This Capability Is Ranked First

Restricting the comparison between NRPyLeakage and a `bns_nurates`-class
provider to major interaction channels that NRPyLeakage **completely lacks**,
inelastic ν−e± (neutrino–electron/positron) scattering is the **most
important** one. The reasons:

- It is entirely absent from the current kernel: the only scattering carried
  is elastic neutral-current scattering on free neutrons and protons. Electrons
  and positrons are not a scattering target at all.
- Among the channels that are completely absent, it has the broadest reach:
  every aggregate species scatters on the thermal lepton background, so its
  absence is felt across the whole `nue`/`anue`/`nux` budget, not just in one
  charged-current family.
- It is not an accuracy correction to an existing channel (see
  [Weak-Magnetism, Recoil, Phase-Space, And Mean-Field Corrections](weak-magnetism-recoil-phase-space-mean-field.md))
  and it is not a situational threshold effect (see
  [Distinct Nucleon Decay Contributions](distinct-nucleon-decay-contributions.md));
  it is a missing interaction class with a spectral character the current grey
  model cannot represent at all.

The other two capability leaves close gaps in channels the current model
already contains. This leaf closes the one gap that is both a missing channel
and a missing interaction class.

## Purpose

This page is a to-do for adding inelastic neutrino-electron and
neutrino-positron scattering to NRPyLeakage, to a level comparable to what a
`bns_nurates`-class rate code offers. It names the physics, the gap it closes
in the current grey leakage model, the architectural and contract changes the
change would force, and the validation that must pass before the channel could
be claimed. Repo-local source, headers, and tests remain the authority for what
exists today.

## Current State (ground truth)

The current leakage scattering channels are elastic neutral-current scattering
on free neutrons and free protons, for all three aggregate species. In the
checked-in kernel the scattering opacities follow the Ruffert et al. (1996)
parametrization, with the prefactor
`NRPyLeakage_N_A * NRPyLeakage_sigma_0 * T^2 * rho_cgs / m_e_c2^2` and the
`alpha`/`C_V` weak-coupling factors, plus neutron and proton blocking factors
of the form `1/(1 + (2/3) fmax(mu_n/T, 0))` and
`1/(1 + (2/3) fmax(mu_p/T, 0))`, and they enter `ghl_neutrino_opacities` slot
`[0]` (number) and slot `[1]` (energy). The checked-in process set is beta
capture, pair annihilation, transverse plasmon decay, and nucleon-nucleon
bremsstrahlung, plus that elastic nucleon scattering.

There is no scattering on electrons or positrons in the current source.
Electrons and positrons appear only as the charged-current partners
(`e^- + p <-> n + nu_e` and the charge-conjugate) and as the source of the
pair and plasmon channels; they are not a scattering target. The electron
chemical potential `mu_e` and temperature `T` are available at the leakage
boundary, but no electron/positron number density, occupation, or thermal
scattering kernel is computed or consumed.

Ground truth for the current boundary:

- [`NRPyLeakage_compute_neutrino_opacities.c`](../../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities.c)
- [`NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c`](../../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c)
- [`ghl_radiation.h`](../../../../GRHayL/include/ghl_radiation.h)
- [Physics And EOS Contract](../physics-and-eos-contract.md)
- [Implementation Flow](../implementation-flow.md)

## The Capability To Investigate

`bns_nurates`-class codes typically carry inelastic neutrino-electron and
neutrino-positron scattering, for example

- `nu_i + e^- <-> nu_i + e^-`,
- `nu_i + e^+ <-> nu_i + e^+`,
- `anti-nu_i + e^+ <-> anti-nu_i + e^+`,
- `anti-nu_i + e^- <-> anti-nu_i + e^-`,

for the three aggregate species, with an additional charged-current (W-exchange)
contribution for the electron flavor on top of the neutral-current channel. The
key physical fact is that these channels are inelastic: the neutrino changes
energy in the electron rest frame, so the process both transfers momentum and
redistributes neutrino energy. In a
grey, energy-averaged model that carries only a total energy density and an
implied mean energy, inelastic scattering has no exact local representation.

## Why This Is Not A Local Leakage Extension

This is the central difficulty and must be stated up front. The current
NRPyLeakage model is grey and energy-averaged: it does not evolve a neutrino
spectrum and it does not carry a neutrino mean energy as an independent
quantity. Inelastic scattering is, by construction, a spectral process: its
effect on the energy budget depends on the shape of the incoming and outgoing
neutrino distributions and on the relative angle between neutrino and target
electron. A grey leakage model cannot represent that effect exactly. Any grey
treatment would be an effective closure: a calibrated energy-transfer term or
an effective opacity that assumes a mean energy and an angular distribution the
model does not actually carry.

This is a general limitation of ordinary grey leakage, not a peculiarity of
this one channel: **ordinary grey leakage cannot self-consistently calculate
processes that depend on the existing non-equilibrium neutrino radiation
field.** The most important examples are

1. **full energy redistribution from inelastic neutrino–electron scattering**
   (this page): the energy transfer depends on the incoming neutrino spectrum
   and the relative neutrino–electron angle, neither of which the grey model
   carries;
2. **inverse pair reactions** (e.g., `nu + anti-nu <-> e^- + e^+`): the
   reaction rate depends on the coexisting neutrino and antineutrino occupation
   distributions, not only on the thermodynamic state;
3. **absorption of neutrinos emitted elsewhere**: a local grey absorption rate
   multiplies the local opacity by a locally inferred radiation intensity, so
   it cannot account for neutrinos whose energy and direction are set by an
   emission region distant from the absorption point.

The reason all three are out of reach of grey leakage is the same: **grey
leakage does not evolve neutrino spectra, densities, and fluxes.** Its local
rates are built from the thermodynamic state (T, density, composition,
degeneracies) and an interpolation between free-streaming and diffusion limits,
so any process whose rate depends on the radiation field itself—its spectrum,
its local density, or its flux from a distant source—can only enter as an
approximate emissivity or opacity, never as a self-consistent calculation.

Consequences for the to-do:

- A grey inelastic-electron channel is scientifically an approximation, not an
  exact leakage term. The page must say what information is discarded, in the
  spirit of the existing "effective rates" and "assumptions and limitations"
  sections of the physics contract.
- If a faithful treatment is required, the change is not a leakage extension at
  all: it requires a mean-energy or multigroup quantity to be transported or
  supplied, which is a larger architecture decision that changes the
  thermodynamic-state-only contract of the current six-file kernel.
- The emissivity-versus-inverse asymmetry also applies: an inelastic scattering
  kernel is harder to fold into an inverse absorption than a pure emissivity
  would be. Do not present the scattering as if it were a local source.

## What The Change Would Need

Investigating this capability would require, at minimum:

1. A decision on the representation: an effective grey closure with an explicit
   stated mean-energy and angular assumption, versus a transport of a grey
   mean energy, versus deferral to a multigroup framework. This decision gates
   everything below.
2. Electron and positron number densities and occupations at the leakage
   boundary. The current EOS callback returns `mu_e`, `mu_n`, `mu_p`, `muhat`,
   `X_n`, and `X_p`; it does not return electron/positron number densities. A
   new channel would need those quantities, either from an extended EOS output
   or from a thermally consistent estimate that is itself documented and
   validated.
3. Cross sections or a table for the inelastic channel, in the units and
   conventions of the leakage boundary, with blocking, thermal, and
   weak-magnetism treatment stated. These belong to
   [Weak-Magnetism, Recoil, Phase-Space, And Mean-Field Corrections](weak-magnetism-recoil-phase-space-mean-field.md)
   once a representation is chosen.
4. A slot or output decision: which `ghl_neutrino_opacities` entries the new
   opacity occupies, whether it is number-conserving (scattering) or
   number-exchanging (charged-current charge exchange), and how it composes
   with the existing elastic nucleon scattering opacity.
5. A source-term decision: inelastic scattering does not, in the isoenergetic
   limit, change number; its energy effect in the current grey model must be
   either represented as an effective energy-exchange term or explicitly
   deferred. The `R_source`/`Q_source` sign conventions in the physics contract
   must not be silently altered.
6. A unit and finite-handling path consistent with the existing
   `EnsureFinite`/`robust_isfinite` usage and the conversion macros owned by
   [`ghl_nrpyleakage.h`](../../../../GRHayL/include/ghl_nrpyleakage.h). Do not
   invent a new conversion set; route to the existing `NRPyLeakage_units_*`
   macros.

## Architectural And Contract Impact

The change would touch, or would have to be audited against, the following:

- `ghl_neutrino_opacities` slot semantics: the two-slot per-species contract is
  currently number/energy. An inelastic channel that changes energy but not
  number must be placed in a way that does not corrupt the diffusion time or
  the luminosity suppression that consume slot `[1]`.
- The EOS boundary: new electron/positron density inputs would either extend
  `ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T` or add a new callback,
  which is a contract change with downstream EOS and GRHayLib consequences.
  See [API And Data](../api-and-data.md) for the current callback and the
  `GHL_DISABLE_HDF5` guard.
- The build manifest
  [`make.code.defn`](../../../../GRHayL/Neutrinos/NRPyLeakage/make.code.defn):
  a new channel is a formula change inside an existing generated block or a new
  generated block. Per the regeneration rule in
  [Generator Provenance](../generator-provenance.md), it should be derived
  from a symbolic expression and compared against the current kernel, not
  hand-written into a `tmp_*` block.
- The public header `ghl_nrpyleakage.h`: any new constant or helper belongs
  here, with the `NRPyLeakage_` prefix, alongside the existing constant family.

## Validation Requirements

This capability is not validated until it can be isolated, bounded, and
routed:

- A standalone test that sets the new channel on and off and shows the current
  six-file kernel results are byte-for-byte unchanged when the channel is off
  (exact-equivalence preservation).
- A small optically thin or single-state case with a hand or externally
  computed inelastic-electron opacity or energy-transfer term, replayed as a
  fixture in the existing `Unit_Tests/unit_test_nrpyleakage_*.c` harness, with
  the comparison result actually checked (the current harness discards
  `ghl_pert_test_fail` results; see the assertion gap in
  [Tests And Fixtures](../tests-and-fixtures.md)).
- A bounded statement of the grey approximation: the assumed mean energy, the
  angular assumption, and the density/temperature regime where the closure is
  expected to be acceptable, plus where it is expected to fail.
- A unit-system check that the new term uses the same
  geometric-to-cgs conversion set as the rest of the leakage code.
- A cross-code diagnostic, not an equivalence claim, against a
  `bns_nurates`-class or independent inelastic-electron reference, with the
  two models' difference attributed to the known grey approximation rather
  than to a sign or unit bug.

## Open Questions And Stop Conditions

- Which representation is acceptable (grey closure versus mean-energy versus
  multigroup)? Until this is answered, the channel should not be coded.
- Does the target application actually need inelastic-electron scattering at
  the densities and temperatures where it matters, or is it deferred physics?
  The existing whitepaper-level audit in this codebase classifies inelastic
  scattering as requiring a spectral/multigroup treatment; that classification
  should be reconfirmed against the intended GRHayL use.
- Is the electron/positron number density available from the EOS, or must it be
  estimated? If estimated, the estimate must be validated independently.
- Stop if the only correct treatment requires a transported mean energy or a
  radiation-field-dependent rate: that is a transport change, not a leakage
  formula, and it is outside the scope of this folder.

## References

Repo-local implementation authority (current behavior):

- [`GRHayL/Neutrinos/NRPyLeakage/`](../../../../GRHayL/Neutrinos/NRPyLeakage/)
- [`GRHayL/include/ghl_radiation.h`](../../../../GRHayL/include/ghl_radiation.h)
- [`GRHayL/include/ghl_nrpyleakage.h`](../../../../GRHayL/include/ghl_nrpyleakage.h)
- [`GRHayL/include/ghl_eos_functions.h`](../../../../GRHayL/include/ghl_eos_functions.h)
- [`Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c`](../../../../Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c)
- [`Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c`](../../../../Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c)

External capability and physics references (recheck before relying on any
detail; the `bns_nurates` name is a lookup seed, not an asserted URL):

- `bns_nurates` binary-neutron-star neutrino-rate code base (capability
  reference for inelastic neutrino-electron and neutrino-positron scattering)
- Standard inelastic neutrino-electron scattering cross-section literature
  (charged- and neutral-current `nu`/`anti-nu` on `e`/`e^+`), for the
  parameterization and the energy-transfer structure.
