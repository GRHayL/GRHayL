# Distinct Nucleon Decay / Capture Channels (Future Work)

Status: future "to do" work. This page is not a record of implemented or
validated functionality. No channel separation beyond the current lumped
nucleon coupling is present in the checked-in NRPyLeakage source, and none is
claimed here.

## Purpose And Position In The Gap

`bns_nurates`-class codes also add **nucleon decay / inverse decay** as
separately parameterized channels relative to the current single beta rate. In
the capability-gap ranking of [The Capability Gap, Stated](index.md#the-capability-gap-stated),
this addition is **generally more situational** than the inelastic ν−e±
scattering that the current model completely lacks: the per-channel structure
is most important near threshold and at high degeneracy, and less important in
the non-degenerate low-temperature limit. It refines channels the current
kernel already contains rather than adding a missing interaction class.

This page is a to-do for representing the charged-current weak channels on free
nucleons as distinct contributions, each with its own phase space, blocking,
weak-magnetism, and sign structure, to a level comparable to a
`bns_nurates`-class rate code. It names the channels, the gap it closes, the
contract changes the change would force, and the validation that must pass
before the change could be claimed.

## Current State (ground truth)

The current leakage charged-current physics is built around two local
approximations and two composition factors:

- Charged-current emission is electron capture on protons
  (`e^- + p -> n + nu_e`) for the `nue` channel and positron capture on
  neutrons (`e^+ + n -> p + anti-nu_e`) for the `anue` channel.
- Charged-current absorption is the charge-conjugate inverse channel for
  `nue` on neutrons and for `anue` on protons.
- The kernel sets `Y_p = Y_e` and `Y_n = 1 - Y_e` for the charged-current and
  scattering factors, and uses the free-nucleon mass fractions `X_n` and `X_p`
  separately for the nucleon-nucleon bremsstrahlung composition factor.
- The beta-process rate in the checked-in kernel is a single per-flavor rate
  built from Fermi-Dirac integral factors (through
  `NRPyLeakage_Fermi_Dirac_integrals`), a shared weak coupling, and a
  blocking factor of the form `1/(1 + (2/3) fmax(mu_n/T, 0))` for neutrons and
  `1/(1 + (2/3) fmax(mu_p/T, 0))` for protons. It carries a single linear
  nucleon blocking factor per flavor; it does not carry a separate
  final-state-lepton blocking, a Fermi-function/Coulomb correction, or a
  per-channel weak-magnetism or recoil term.

The current source therefore distinguishes the neutron and proton blocking
factors, but it does not carry the forward capture and the inverse absorption
as separately parameterized channels with their own thresholds, their own
weak-magnetism coefficients, and their own final-state phase-space factors.
The free nucleon beta-decay channel `n -> p + e^- + anti-nu_e` and its
conjugate `p -> n + e^+ + nu_e` are not an independent rate in the current
kernel; the relevant physics enters through the capture rates and the
equilibrium neutrino degeneracy.

Ground truth for the current boundary:

- [`NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c`](../../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c)
- [`NRPyLeakage_compute_neutrino_luminosities.c`](../../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_luminosities.c)
- [`ghl_nrpyleakage.h`](../../../../GRHayL/include/ghl_nrpyleakage.h)
- [Physics And EOS Contract](../physics-and-eos-contract.md)
- [Implementation Flow](../implementation-flow.md)

## The Capability To Investigate

A `bns_nurates`-class treatment keeps the charged-current weak channels
separate so that each can carry its own physics:

- Electron capture: `e^- + p -> n + nu_e` (raises `Y_e`-relevant neutron
  production; emits `nu_e`).
- Positron capture: `e^+ + n -> p + anti-nu_e`.
- Inverse beta / neutrino absorption: `nu_e + n -> p + e^-`.
- Antineutrino absorption: `anti-nu_e + p -> n + e^+`.
- Free nucleon beta decay: `n -> p + e^- + anti-nu_e` and the conjugate
  `p -> n + e^+ + nu_e`.

Each of these has a distinct phase-space endpoint (set by the neutron-proton
mass gap and the lepton energy), a distinct final-state particle set (which
determines which blocking factor applies), a distinct weak-magnetism
contribution, and a distinct sign in the lepton-number source `R_source`.
Representing them as separate contributions lets the model:

- apply the correct final-state blocking per channel (neutron blocking for
  channels that emit a neutron, proton blocking for channels that emit a
  proton), rather than a single shared factor;
- carry the neutron-proton mass gap `Q_npmass` in the channel thresholds where
  the current kernel consumes `muhat` directly;
- separate the forward and inverse detailed-balance factors explicitly, which
  the grey equilibrium interpolation currently folds into a single
  equilibrium-degeneracy interpolation.

## What The Change Would Need

Investigating this capability would require, at minimum:

1. A per-channel rate table or per-channel generated formulas, one for each of
   the six channels above, each expressed at the leakage boundary in MeV and
   cgs, in the units of the existing kernel.
2. A decision on the neutron-proton mass gap. The current kernel consumes
   `muhat` directly in the equilibrium combination `(mu_e - muhat)/T` and does
   not subtract a rest-mass gap `Q_npmass` in the leakage routine. The
   constants `NRPyLeakage_Q_npmass` and `NRPyLeakage_ZL_Q_npmass` are declared
   in [`ghl_nrpyleakage.h`](../../../../GRHayL/include/ghl_nrpyleakage.h) but
   are not consumed by any of the six built files. Introducing `Q_npmass` into
   channel thresholds is a convention change that must be made deliberately and
   revalidated; see the chemical-potential convention hazard in
   [Physics And EOS Contract](../physics-and-eos-contract.md).
3. Distinct blocking factors per channel, with the final-state particle named
   for each, not a single shared `fmax(mu_N/T, 0)` term.
4. A distinct weak-magnetism coefficient per channel, which is shared with
   [Weak-Magnetism, Recoil, Phase-Space, And Mean-Field Corrections](weak-magnetism-recoil-phase-space-mean-field.md).
5. A lepton-number source assembly that preserves the current
   `R_source = R_{anue}^{eff} - R_{nue}^{eff}` sign convention and the
   electron-fraction meaning, while allowing each channel to contribute its own
   signed rate. The `nux` factor of four and the absence of a heavy-lepton
   number rate must be preserved unless deliberately changed.
6. A detailed-balance or equilibrium-consistency statement for each
   forward/inverse pair, so that the grey equilibrium interpolation remains the
   fixed point of the assembled channels.

## Why This Is More Than A Constant Change

The tempting shortcut is to treat this as adding a `Q_npmass` term and a
couple of prefactors. It is not. The current kernel's beta rate is a single
symbolic expression whose `tmp_*` graph assumes a particular channel structure.
Separating the channels changes the structure of that expression, not just its
constants. Per the regeneration rule in
[Generator Provenance](../generator-provenance.md), the per-channel forms
should be derived from symbolic expressions and compared against the current
kernel, then recombined. Hand-editing the `tmp_*` block to split channels
destroys the provenance link to the ancestral generator and makes the
exact-equivalence property of the current kernel impossible to check.

## Architectural And Contract Impact

The change would touch, or would have to be audited against, the following:

- The charged-current rate block inside
  [`NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c`](../../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c)
  and the matching luminosity block in
  [`NRPyLeakage_compute_neutrino_luminosities.c`](../../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_luminosities.c).
- The `R_source`/`Q_source` assembly and the heavy-lepton factor of four in the
  combined source routine.
- The constants in [`ghl_nrpyleakage.h`](../../../../GRHayL/include/ghl_nrpyleakage.h):
  per-channel weak-magnetism and phase-space coefficients, and a decision on
  whether `Q_npmass` becomes live.
- The EOS boundary if a channel needs a composition quantity not returned by
  `ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T` (for example, separate
  neutron and proton number fractions distinct from `X_n`/`X_p` and
  `Y_n`/`Y_p`). Route that through
  [API And Data](../api-and-data.md).
- The build manifest
  [`make.code.defn`](../../../../GRHayL/Neutrinos/NRPyLeakage/make.code.defn)
  and the unit tests and CI route in
  [Tests And Fixtures](../tests-and-fixtures.md).

## Validation Requirements

This capability is not validated until:

- A per-channel on/off test shows each channel can be isolated, and that the
  current combined result is reproduced when all channels are enabled in the
  current-limiting regime (low temperature, non-degenerate nucleons) so that
  the new structure reduces to the old one.
- A detailed-balance test verifies each forward/inverse pair at equilibrium,
  using the existing equilibrium-degeneracy definition, so that the assembled
  channels do not move the weak-equilibrium fixed point.
- A sign and lepton-number test verifies that the assembled `R_source`
  preserves the current electron-fraction meaning and the heavy-lepton
  factor-of-four behavior.
- A blocking test checks the neutron and proton blocking factors against an
  independent calculation at high degeneracy.
- A unit-system and finite-handling check consistent with the existing
  `EnsureFinite`/`robust_isfinite` usage and the
  [`ghl_nrpyleakage.h`](../../../../GRHayL/include/ghl_nrpyleakage.h)
  conversion macros.
- A cross-code diagnostic, not an equivalence claim, against a
  `bns_nurates`-class or independent per-channel reference, with the
  difference attributed to the stated grey and blocking approximations.

## Open Questions And Stop Conditions

- Does the target application need per-channel resolution, or is the current
  single beta rate adequate for the intended density and temperature range?
  Per-channel structure is most important near threshold and at high
  degeneracy; it is less important in the non-degenerate low-temperature
  limit.
- Is `Q_npmass` to be introduced as a live threshold, or is the current
  `muhat`-direct convention to be preserved? This is a convention decision
  that interacts with the EOS chemical-potential contract; it must not be made
  silently.
- Which channels are forward-emitting and which are inverse-absorbing in the
  grey relaxation model, and how does each map to an opacity slot?
- Stop and route to the architecture leaf if the per-channel structure forces a
  new EOS output or a new radiation struct field rather than a rate-formula
  change.

## References

Repo-local implementation authority (current behavior):

- [`GRHayL/Neutrinos/NRPyLeakage/`](../../../../GRHayL/Neutrinos/NRPyLeakage/)
- [`GRHayL/include/ghl_nrpyleakage.h`](../../../../GRHayL/include/ghl_nrpyleakage.h)
- [`GRHayL/include/ghl_radiation.h`](../../../../GRHayL/include/ghl_radiation.h)
- [`GRHayL/include/ghl_eos_functions.h`](../../../../GRHayL/include/ghl_eos_functions.h)
- [`Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c`](../../../../Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c)

External capability and physics references (recheck before relying on any
detail; the `bns_nurates` name is a lookup seed, not an asserted URL):

- `bns_nurates` binary-neutron-star neutrino-rate code base (capability
  reference for per-channel charged-current and free-decay treatment).
- Standard charged-current weak-interaction rate literature for the per-channel
  phase space, the neutron-proton mass gap, and the detailed-balance
  relations.
