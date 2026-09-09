# NRPyLeakage Future-Work Wiki: Toward `bns_nurates`-Class Capabilities

## Purpose

This is a future-work wiki folder, not a record of implemented or validated
functionality. It catalogs distinct, separately scoped enhancements that
NRPyLeakage could investigate so its rate, opacity, and source-term formulas
reach a higher fidelity than the grey leakage model it currently implements.
The capabilities are named after what an established binary-neutron-star
neutrino-rate code base such as `bns_nurates` offers. `bns_nurates` is used as
a capability reference only; it is not assumed to be a repo-local dependency
and it is not the source of the current GRHayL leakage formulas. None of the
work described here is complete.

## The Capability Gap, Stated

Compared with a `bns_nurates`-class provider, the gap in the current
NRPyLeakage surface breaks down as follows:

- **Restricting the comparison to major interaction channels that NRPyLeakage
  completely lacks, inelastic ν−e± (neutrino–electron/positron) scattering is
  the most important one.** The current kernel carries elastic
  neutral-current scattering on free neutrons and protons only; there is no
  scattering on electrons or positrons at all. See
  [Inelastic Neutrino-Electron/Positron Scattering](inelastic-neutrino-electron-scattering.md).
- **`bns_nurates` also adds nucleon decay / inverse decay** as separately
  parameterized charged-current channels (forward capture, inverse absorption,
  and free beta decay, each with its own phase space, blocking, and
  weak-magnetism factors). **That is generally more situational:** it is most
  important near threshold and at high degeneracy, and less important in the
  non-degenerate low-temperature limit. See
  [Distinct Nucleon Decay Contributions](distinct-nucleon-decay-contributions.md).
- **The larger advantage of `bns_nurates` is that it provides substantially
  more accurate physics for several of the dominant interactions that both
  providers already contain:** weak-magnetism, recoil, phase-space
  (blocking/Fermi-function), and relativistic mean-field corrections to the
  shared scattering and charged-current channels. See
  [Weak-Magnetism, Recoil, Phase-Space, And Mean-Field Corrections](weak-magnetism-recoil-phase-space-mean-field.md).
- **The transverse plasmon decay that NRPyLeakage includes is not one of the
  dominant processes in BNS mergers.** It is one of the checked-in thermal
  channels, but it does not set the BNS-neutrino budget, and its absence from a
  rate-code comparison is not what makes `bns_nurates` more capable.
- **No interaction is physically valid only for M1.** Leakage can include
  approximate emissivities or opacities for all of them; the grey leakage
  model is not restricted to a subset of the interaction list for
  physics-validity reasons.
- **However, ordinary grey leakage cannot self-consistently calculate
  processes that depend on the existing non-equilibrium neutrino radiation
  field.** In particular: full energy redistribution from inelastic
  neutrino–electron scattering, inverse pair reactions, and absorption of
  neutrinos emitted elsewhere. Grey leakage does not evolve neutrino spectra,
  densities, and fluxes, so its local rates are built from thermodynamic state
  and an equilibrium-interpolated degeneracy, not from the radiation field
  that those processes act on.

## Current Model (ground truth)

The current leakage model is grey and local: it computes free emission rates,
an energy-averaged opacity, a path-of-least-resistance optical depth, and an
interpolation between free-streaming and diffusion-limited rates. Its
checked-in processes are electron and positron capture on free nucleons,
electron-positron pair annihilation, transverse plasmon decay, and
nucleon-nucleon bremsstrahlung, plus absorption and elastic scattering on
free neutrons and protons.

### Reading `bns_nurates` as a capability reference

`bns_nurates` is treated the way the Neutrinos KB treats the ancestral
`Tabulated_EOS_IllinoisGRMHD` repository: as provenance and a capability
reference, not as the authority for current GRHayL signatures or behavior. Its
role in this folder is to name the target microphysics named in
[The Capability Gap, Stated](#the-capability-gap-stated).

Before any incorporation work, recheck the actual `bns_nurates` repository,
its process list, its cross-section parameterizations, its units and
conventions, and its test or verification route. Do not assume the names used
in this folder match its internal symbols, and do not assume its numbers are
drop-in. No stable repository URL is asserted here; treat the `bns_nurates`
name as a lookup seed and confirm the authoritative source before relying on
any external detail.

### What "consistent with NRPyLeakage and GRHayL" means

An incorporation is GRHayL-consistent only if it preserves, or deliberately and
documentedly changes, each of the following and then revalidates them:

- the three-species `nue`/`anue`/`nux` ordering and the `[0]`/`[1]`
  number/energy slot contract of `ghl_neutrino_opacities`;
- the `ghl_error_codes_t` return convention and the before-writeback failure
  boundary;
- the `ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T` EOS callback and
  the `GHL_DISABLE_HDF5` guard on the three EOS-dependent routines;
- the geometric-to-cgs unit conversions owned by
  [`ghl_nrpyleakage.h`](../../../../GRHayL/include/ghl_nrpyleakage.h);
- the `NRPyLeakage_*` public spelling and the radiation-struct `ghl_*`
  spelling; and
- the build manifest and the unit-test and CI route.

## Folder Map

Each leaf below is a separate "to do" workstream. None is implemented.

| Leaf | Future capability | Primary risk if done carelessly |
| --- | --- | --- |
| [Inelastic neutrino-electron/positron scattering](inelastic-neutrino-electron-scattering.md) | Add energy-exchanging `nu_i`/`anti-nu_i` scattering on thermal electrons and positrons. | Implies a spectral/mean-energy structure that grey leakage does not carry; risks inventing a fake grey channel. |
| [Distinct nucleon decay contributions](distinct-nucleon-decay-contributions.md) | Separate charged-current capture and decay channels on free neutrons from those on free protons, with distinct blocking, phase space, and weak-magnetism factors. | Lumping neutrons and protons hides channel-specific thresholds, blocking, and sign structure. |
| [Weak-magnetism, recoil, phase-space, and mean-field corrections](weak-magnetism-recoil-phase-space-mean-field.md) | Upgrade the shared scattering and charged-current channels beyond the current leading-order approximation. | A correction folded into a grey rate without a stated averaging rule is unquantifiable and untestable. |

## Cross-Cutting Ground Truth

- Source: `GRHayL/Neutrinos/NRPyLeakage/` (the six built files named in
  [Implementation Flow](../implementation-flow.md)).
- Public radiation structs:
  [`GRHayL/include/ghl_radiation.h`](../../../../GRHayL/include/ghl_radiation.h)
- Public leakage declarations, constants, and unit conversions:
  [`GRHayL/include/ghl_nrpyleakage.h`](../../../../GRHayL/include/ghl_nrpyleakage.h)
- Error codes: [`GRHayL/include/ghl.h`](../../../../GRHayL/include/ghl.h)
- EOS callback:
  [`GRHayL/include/ghl_eos_functions.h`](../../../../GRHayL/include/ghl_eos_functions.h)
  and
  [`GRHayL/include/ghl_eos_functions_declaration.h`](../../../../GRHayL/include/ghl_eos_functions_declaration.h)
- Build manifest:
  [`GRHayL/Neutrinos/NRPyLeakage/make.code.defn`](../../../../GRHayL/Neutrinos/NRPyLeakage/make.code.defn)
- Tests: `Unit_Tests/nrpyleakage_main.h` and
  `Unit_Tests/unit_test_nrpyleakage_*.c`
- Ancestral derivation and generator evidence:
  [Generator Provenance](../generator-provenance.md)

## Shared Rules For This Folder

- Frame every capability as future work. Use "would", "should", and "to do";
  do not use "adds", "provides", or "validates" to describe current behavior.
- Keep physics statements bounded to what the current source and the cited
  derivation evidence support. A correction that is physically real is not a
  claim that it is present in the current leakage code.
- Every capability leaf must state: the gap it closes, the architectural and
  contract changes it forces, the validation it requires, and the specific
  current-surface files and tests it would touch.
- Treat any cross-check against `bns_nurates` as a cross-code diagnostic, not
  an equivalence claim: NRPyLeakage and a `bns_nurates`-class provider are
  different physical models unless a common backend is built.
- Do not copy `bns_nurates` wholesale. Incorporation is per-capability, per the
  consistency rules above, with the current exact-equivalence kernel preserved
  as the fallback.
