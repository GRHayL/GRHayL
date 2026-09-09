# Weak-Magnetism, Recoil, Phase-Space, And Mean-Field Corrections (Future Work)

Status: future "to do" work. This page is not a record of implemented or
validated functionality. None of the corrections named here is present in the
checked-in NRPyLeakage source, and none is claimed here. The current kernel
uses a leading-order approximation for the shared scattering and
charged-current channels; the corrections below are candidates to investigate.

## Purpose And Position In The Gap

In the capability-gap ranking of [The Capability Gap, Stated](index.md#the-capability-gap-stated),
this leaf is where **the larger advantage of `bns_nurates` lives**: it
provides **substantially more accurate physics for several of the dominant
interactions that both providers already contain**. Unlike the inelastic
ν−e± scattering that NRPyLeakage completely lacks, and unlike the more
situational nucleon-decay channel separation, the weak-magnetism, recoil,
phase-space (blocking/Fermi-function), and relativistic mean-field corrections
in this leaf upgrade channels that are already in both models—elastic
neutrino–nucleon scattering and the charged-current beta channels—which is
precisely where the dominant BNS-neutrino interactions sit. The current
kernel carries these channels only at leading order; the corrections below are
what a `bns_nurates`-class rate code adds on top.

This page is a to-do for upgrading the shared neutrino-nucleon and
neutrino-lepton channels beyond the current leading-order approximation, to a
level comparable to a `bns_nurates`-class rate code. It separates four
distinct correction families: weak magnetism, recoil, phase-space
(blocking and Fermi-function) effects, and relativistic mean-field effects. For
each, it names what the current kernel does, what the correction adds, the
contract changes the change would force, and the validation that must pass
before the correction could be claimed.

## Current State (ground truth)

The current kernel is a leading-order, grey, energy-averaged approximation.
Relevant facts about the shared channels:

- The scattering channels are elastic neutral-current on free neutrons and
  protons, following the Ruffert et al. (1996) parametrization, with the
  Thomson-like prefactor
  `NRPyLeakage_N_A * NRPyLeakage_sigma_0 * T^2 * rho_cgs / m_e_c2^2` and
  weak-coupling factors built from `NRPyLeakage_alpha`,
  `NRPyLeakage_C_V`, `NRPyLeakage_gamma_0`, and `NRPyLeakage_alpha_fs`, plus
  neutron and proton blocking factors `1/(1 + (2/3) fmax(mu_n/T, 0))` and
  `1/(1 + (2/3) fmax(mu_p/T, 0))`. The checked-in kernel does not expose a
  separate, independently tuned weak-magnetism term on the nucleon nor a
  nucleon-recoil factor; any weak-magnetism content is only what the fixed
  `C_V`-based neutral-current coefficient already encodes.
- The charged-current beta rate is a single symbolic expression with a
  blocking factor `1/(1 + (2/3) fmax(mu_n/T, 0))` (neutron) or
  `1/(1 + (2/3) fmax(mu_p/T, 0))` (proton), a shared Fermi factor, and the
  `muhat`-direct equilibrium combination `(mu_e - muhat)/T`. The checked-in
  charged-current block does not expose a separate weak-magnetism term, a
  recoil factor, or a relativistic mean-field shift; it uses the leading-order
  Fermi coupling and a single linear blocking factor.
- The pair and plasmon channels and the nucleon-nucleon bremsstrahlung channel
  are the leading-order thermal-pair and bremsstrahlung forms; they are not
  the target of this page except where a correction touches the electron
  chemical potential or the density normalization.

The relevant constants live in
[`ghl_nrpyleakage.h`](../../../../GRHayL/include/ghl_nrpyleakage.h); the
relevant rate blocks live in the six built files named in
[Implementation Flow](../implementation-flow.md).

## The Four Correction Families

### Weak Magnetism

Weak magnetism is the nucleon magnetic-moment contribution to the
charged-current (and, for some conventions, to the neutral-current) weak
current. For neutron targets it is the dominant correction to the charged-
current rate and is parameterized by the neutron magnetic moment, commonly
written with a factor proportional to `mu_n` (in nuclear magnetons) relative to
the Fermi coupling. The current kernel carries no explicit weak-magnetism term
on the charged-current channel.

To investigate:

- Identify the exact weak-magnetism factor for each charged-current channel
  (capture, absorption, and free decay) in the convention of the leakage
  boundary, including the sign and the energy dependence.
- Determine whether the same or a different weak-magnetism factor applies to
  the elastic scattering channel on neutrons versus protons, and whether the
  current `C_V`-based neutral-current form already absorbs part of it.
- Decide the averaging: in a grey model the weak-magnetism contribution is
  folded into an energy-averaged rate. State the mean energy or the
  averaging weight used, because the factor is energy dependent.
- The neutron-proton mass gap and the weak-magnetism coefficient both touch the
  neutron-proton chemical-potential combination; this interacts with the
  chemical-potential convention hazard in
  [Physics And EOS Contract](../physics-and-eos-contract.md) and with the
  per-channel structure in
  [Distinct Nucleon Decay Contributions](distinct-nucleon-decay-contributions.md).

### Recoil

Recoil is the correction from the finite nucleon mass in the two-body
kinematics: the final nucleon does not remain at rest, so the available phase
space and the energy transfer differ from the infinite-mass limit. For
charged-current and elastic scattering channels this is a `1/m_N`-suppressed
correction that becomes more important at higher neutrino energy and for
lighter targets. The current kernel uses the infinite-mass, leading-order
form.

To investigate:

- Determine the recoil factor for each channel in the convention of the
  leakage boundary, and the range of neutrino energy over which it is
  non-negligible at the densities and temperatures of interest.
- Decide the averaging, as with weak magnetism: the recoil factor is energy
  dependent and must be folded into the grey rate with a stated mean energy.
- Confirm whether recoil affects the number opacity, the energy opacity, or
  both, and whether it shifts the equilibrium target used by the grey
  relaxation.

### Phase Space (Blocking And Fermi Functions)

Phase-space effects here are the final-state occupancy corrections: Pauli
blocking of the final nucleon and lepton, and the Fermi-function (Coulomb)
correction to the emitted or absorbed lepton. The current kernel carries a
linear blocking factor `1/(1 + (2/3) fmax(mu_N/T, 0))` on the charged-current
channel and uses Fermi-Dirac integrals (keys 0-5) for the thermal
populations. The Fermi-function/Coulomb correction and a more complete
final-state-lepton blocking are not present in the checked-in charged-current
block.

To investigate:

- Determine the full Fermi-function factor for the emitted/absorbed electron or
  positron, including the sign for attraction versus repulsion, and the
  regime where the approximation used here (or a tabulated Fermi function) is
  valid.
- Determine the complete final-state-lepton blocking in addition to the
  nucleon blocking, for the channels where the lepton is the final state.
- Decide whether the existing Fermi-Dirac integral helper
  `NRPyLeakage_Fermi_Dirac_integrals` (keys 0-5) is sufficient for the new
  occupancy factors or whether additional moments are required, and route any
  new key through the invalid-key error contract in
  [Implementation Flow](../implementation-flow.md).

### Relativistic Mean-Field Effects

Relativistic mean-field effects are the modifications to the nucleon and
lepton energies, effective masses, and densities from the dense-matter mean
field: vector and scalar potentials that shift the nucleon chemical
potentials, modify the neutron-proton asymmetry potential, and change the
effective nucleon mass that enters the phase space and the cross-section
normalization. The current kernel uses the bare `mu_n`, `mu_p`, `muhat`,
`X_n`, and `X_p` from the EOS callback and does not apply an explicit
mean-field correction to the rates.

To investigate:

- Determine which mean-field quantities the EOS already provides (some tabulated
  EOS return effective masses or potentials) and which would have to be added
  to the callback. Route any new EOS output through
  [API And Data](../api-and-data.md).
- Decide whether the mean-field correction is folded into the rate as an
  effective-mass and potential shift, or whether it changes the density and
  chemical-potential normalization upstream of the rate.
- State the density regime where the mean-field correction is expected to
  matter (near and above nuclear density) and where it is expected to vanish.
- Confirm that the mean-field shift is consistent with the EOS convention for
  `muhat`; a mean-field correction that moves `muhat` independently of the EOS
  breaks the chemical-potential contract.

## What All Four Corrections Share

- Each is energy or density dependent and therefore must be averaged for a
  grey model. The averaging rule (mean energy, density weight, or a tabulated
  mean) must be stated per correction; otherwise the correction is not
  reproducible.
- Each belongs to the symbolic kernel, not to the C `tmp_*` block directly. Per
  the regeneration rule in [Generator Provenance](../generator-provenance.md),
  the corrected forms should be derived symbolically, generated into a
  disposable directory, compared against the current kernel, and recombined
  with the GRHayL ABI, EOS, error, unit, and finite-handling contracts reapplied.
- Each preserves the current exact-equivalence property: with the correction
  disabled (coefficient set to its leading-order value or the correction term
  removed), the current six-file kernel results must be reproduced unchanged.

## Architectural And Contract Impact

The change would touch, or would have to be audited against, the following:

- The scattering and charged-current rate blocks inside
  [`NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c`](../../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c)
  and the matching luminosity block in
  [`NRPyLeakage_compute_neutrino_luminosities.c`](../../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_luminosities.c).
- The opacity blocks in
  [`NRPyLeakage_compute_neutrino_opacities.c`](../../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities.c)
  for the scattering corrections.
- The constants in [`ghl_nrpyleakage.h`](../../../../GRHayL/include/ghl_nrpyleakage.h):
  weak-magnetism coefficients, recoil factors, Fermi-function parameters, and
  mean-field potential inputs, each with the `NRPyLeakage_` prefix.
- The EOS boundary if a mean-field or effective-mass quantity is required, via
  `ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T` or a new callback.
- The Fermi-Dirac integral helper and its key/error contract if new occupancy
  moments are required.
- The build manifest and the unit tests and CI route.

## Validation Requirements

A correction is not validated until it can be isolated, bounded, and routed:

- A per-correction on/off test that shows the current kernel is reproduced
  unchanged when the correction is in its leading-order (disabled) form.
- A bounded statement of the averaging rule for each correction and the
  density and temperature regime where the correction is expected to be
  acceptable, plus where it is expected to be negligible or invalid.
- A regime test that confirms each correction is negligible in the
  low-density, low-temperature, non-degenerate limit, so that the corrected
  kernel reduces to the current one outside its intended regime.
- A detailed-balance and equilibrium-consistency check that the correction does
  not move the weak-equilibrium fixed point.
- A unit-system and finite-handling check consistent with the existing
  `EnsureFinite`/`robust_isfinite` usage and the
  [`ghl_nrpyleakage.h`](../../../../GRHayL/include/ghl_nrpyleakage.h)
  conversion macros.
- A cross-code diagnostic, not an equivalence claim, against a
  `bns_nurates`-class or independent corrected reference, with the difference
  attributed to the stated grey and averaging approximations.

## Open Questions And Stop Conditions

- Which corrections are material for the intended application and regime?
  Weak magnetism is usually the largest charged-current correction; recoil is
  usually smallest; mean-field matters only near and above nuclear density.
  A correction that is negligible in the target regime should be deferred.
- Does the EOS provide the mean-field and effective-mass quantities, or must
  they be estimated? If estimated, the estimate must be validated independently
  and routed through the EOS contract.
- Can the averaging be done with the existing Fermi-Dirac integral helper, or
  does a new moment or a tabulated mean energy have to be introduced? A
  tabulated mean energy is a larger architecture change routed to
  [The Capability Gap, Stated](index.md#the-capability-gap-stated).
- Stop and route to the architecture leaf if a correction requires a quantity
  the current EOS callback does not return or a radiation struct field that
  does not exist.

## References

Repo-local implementation authority (current behavior):

- [`GRHayL/Neutrinos/NRPyLeakage/`](../../../../GRHayL/Neutrinos/NRPyLeakage/)
- [`GRHayL/include/ghl_nrpyleakage.h`](../../../../GRHayL/include/ghl_nrpyleakage.h)
- [`GRHayL/include/ghl_radiation.h`](../../../../GRHayL/include/ghl_radiation.h)
- [`GRHayL/include/ghl_eos_functions.h`](../../../../GRHayL/include/ghl_eos_functions.h)
- [`Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c`](../../../../Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c)
- [`Unit_Tests/unit_test_code_error.c`](../../../../Unit_Tests/unit_test_code_error.c)

External capability and physics references (recheck before relying on any
detail; the `bns_nurates` name is a lookup seed, not an asserted URL):

- `bns_nurates` binary-neutron-star neutrino-rate code base (capability
  reference for the corrected channel set).
- Standard weak-interaction correction literature for weak magnetism, recoil,
  Fermi-function, and relativistic mean-field effects on neutrino-nucleon and
  neutrino-lepton rates.
