# Neutrinos EOS Inputs And Chemical Potentials

> Six-output tabulated EOS callback arguments, MeV chemical-potential conventions, equilibrium neutrino degeneracy, `munu` vs. `mu_e - muhat`, and gradual-underflow requirements. · Status: confirmed
> Up: [Neutrinos](index.md)

## Summary

EOS-dependent NRPyLeakage routines obtain thermodynamic chemistry through `ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T`, receiving `(muhat, mu_e, mu_p, mu_n, X_n, X_p)` in that exact output-pointer order. Leakage uses `mu_e` (including electron rest mass) and `muhat` (`mu_n - mu_p`) in MeV to form equilibrium degeneracies and the paired reaction shift, while reconstructing nucleon blocking exclusively from free-nucleon fractions `X_n` and `X_p` rather than raw `mu_n` or `mu_p`.

## Detail

### Six-Output Tabulated EOS Callback Arguments

After `(eos, rho, Y_e, T)`, `ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T` writes six outputs in exact declaration order (`muhat`, `mu_e`, `mu_p`, `mu_n`, `X_n`, `X_p`):

| Output | Physical meaning | Direct consumption in current leakage source |
| --- | --- | --- |
| `muhat` | Neutron-proton chemical-potential difference (`mu_n - mu_p` in GRHayL's StellarCollapse adapter) | `(mu_e - muhat) / T` in grey equilibrium neutrino degeneracies and reaction shift `q` |
| `mu_e` | Electron chemical potential including electron rest-mass energy (O'Connor-Ott/StellarCollapse convention) | Electron degeneracy `mu_e / T`, equilibrium neutrino degeneracy, beta/pair/plasmon rates |
| `mu_p` | Proton chemical potential | Returned for six-output ABI compatibility; not consumed directly by leakage rates |
| `mu_n` | Neutron chemical potential | Returned for six-output ABI compatibility; not consumed directly by leakage rates |
| `X_n` | Free-neutron mass fraction | Neutron scattering population, transition overlap, bremsstrahlung composition factor |
| `X_p` | Free-proton mass fraction | Proton scattering population, transition overlap, bremsstrahlung composition factor |

Chemical potentials and `T` must share MeV units because leakage forms `mu_e / T` and `(mu_e - muhat) / T` without conversion; `Y_e`, `X_n`, and `X_p` are dimensionless. `X_n` and `X_p` are strictly free-nucleon mass fractions, never total nucleon fractions `1 - Y_e` and `Y_e`, so bound nucleons in nuclei are not counted as free scattering or capture targets. GRHayL reads and interpolates StellarCollapse `mu_e` without a zero-point shift; any alternative EOS supplying a rest-mass-subtracted electron chemical potential must add electron rest mass before passing `mu_e` to leakage.

Evidence (scientific/interface requirement): `ghl_eos_functions_declaration.h::ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T`; corroboration: `NRPyEOS_muhat_mue_mup_mun_Xn_and_Xp_from_rho_Ye_T.c`.

### Equilibrium Degeneracies, Stored `munu`, And Unused Mass-Gap Constants

Current leakage source computes equilibrium neutrino degeneracies directly from `mu_e` and `muhat`:

$$\eta_{\nu_e}^{\rm eq} = \frac{\mu_e - \widehat\mu}{T}, \qquad \eta_{\bar\nu_e}^{\rm eq} = -\eta_{\nu_e}^{\rm eq}.$$

Both GRHayL's StellarCollapse table reader and the CompOSE converter define `muhat = mu_n - mu_p`, preserving the neutron-proton rest-mass difference. StellarCollapse documentation notes that `munu = mu_e - mu_n + mu_p = mu_e - muhat`, while warning that some Lattimer-Swesty (LS) tables store an erroneously shifted `munu` dataset. Current NRPyLeakage never reads stored `munu`; it forms `mu_e - muhat` from the six-output callback. By contrast, `GRHayL/EOS/Tabulated/NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T.c` does read stored `munu` to locate beta-equilibrium zero crossings, so the upstream LS table caveat still applies there. Leakage adapters should supply the six callback quantities rather than substituting table `munu`.

`GRHayL/include/ghl_nrpyleakage.h` declares `NRPyLeakage_Q_npmass` and `NRPyLeakage_ZL_Q_npmass`, neither of which is referenced by any current leakage `.c` file. These constants are unrelated to `Q_source`. Do not add or subtract a neutron-proton mass gap $Q$ inside `(mu_e - muhat) / T` or inside nucleon blocking without a paired reaction-model change, and do not reconstruct `muhat` from differently zeroed `mu_n` and `mu_p`. In particular, `muhat(Y_e = 0.5)` is not a rest-mass measurement because even a classical unequal-mass gas has a temperature-dependent kinetic contribution.

Evidence (descriptive behavior): `NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T.c::NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T`; corroboration: `ghl_nrpyleakage.h::NRPyLeakage_Q_npmass`.

### Reaction-Energy Shift, Gradual Underflow, And Adapter Checks

Nucleon blocking and transition overlaps do not consume raw `mu_n` or `mu_p`; instead, `NRPyLeakage_nucleon_blocking.h` inverts free populations `X_n` and `X_p` with a common bare nucleon mass to obtain kinetic degeneracies $\eta_n, \eta_p$, evaluates close-population subtractions using Sterbenz's exact-subtraction theorem, and supplies the shared reaction-energy shift $q = \widehat\mu - T(\eta_n - \eta_p)$ to paired charged-current emission and absorption moments. See [Nucleon Blocking And EOS Conventions](nucleon-blocking-and-eos-conventions.md) for the complete blocking derivation, detailed-balance pairing, producer-energy conventions, and BNS_NURATES/DD2 dense-matter qualification, and [StellarCollapse Table Adapter](../eos/stellarcollapse-table-adapter.md#producer-chemical-potential-energy-references) for legacy SLy4 vs. CompOSE SRO-141 zero points.

No physical clamp is applied to $q$ or $q/T$ because the current EOS callback supplies no mean-field bound; large finite $|q/T|$ remains a conditioning limit. Normalized absorption can divide two representable subnormal Fermi moments, so supported builds require IEEE gradual underflow (`configure` enforces `-no-ftz` on Intel LLVM and probes subnormal behavior); flushing either operand to zero alters a finite physical ratio before final-output sanitization can detect it.

Adapter validation must check all six callback outputs in declaration order, `mu_e / T`, `(mu_e - muhat) / T`, and the derived blocking populations against known thermodynamic states before comparing final rates or opacities. Matching pressure and internal energy alone does not validate the leakage EOS boundary, and matching raw `mu_n/T` or `mu_p/T` no longer tests active blocking.

Evidence (descriptive behavior): `NRPyLeakage_nucleon_blocking.h::nrpyl_compute_population`; corroboration: `NRPyEOS_stellarcollapse.h::NRPyEOS_stellarcollapse_read_table`.

## Sources

- [GRHayL/include/ghl_eos_functions_declaration.h](../../../GRHayL/include/ghl_eos_functions_declaration.h) - `ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T`
- [GRHayL/include/ghl_eos_functions.h](../../../GRHayL/include/ghl_eos_functions.h) - tabulated EOS dispatch declarations
- [GRHayL/include/ghl_nrpyleakage.h](../../../GRHayL/include/ghl_nrpyleakage.h) - `NRPyLeakage_Q_npmass`, `NRPyLeakage_ZL_Q_npmass`
- [GRHayL/EOS/Tabulated/interpolators/NRPyEOS_muhat_mue_mup_mun_Xn_and_Xp_from_rho_Ye_T.c](../../../GRHayL/EOS/Tabulated/interpolators/NRPyEOS_muhat_mue_mup_mun_Xn_and_Xp_from_rho_Ye_T.c) - six-field StellarCollapse interpolator
- [GRHayL/EOS/Tabulated/NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T.c](../../../GRHayL/EOS/Tabulated/NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T.c) - `munu` zero-crossing lookup
- [GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.h](../../../GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.h) - StellarCollapse table reader declaration
- [GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_nucleon_blocking.h](../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_nucleon_blocking.h) - `NRPyLeakage_compute_nucleon_blocking`, Sterbenz subtraction, shifted moments
- [StellarCollapse EOS Table Documentation](https://stellarcollapse.org/equationofstate.html) - `mu_e` rest-mass convention and LS `munu` warning
- [O2scl eos_sn_oo Reference](https://awsteiner.org/code/o2scl/html/class/eos_sn_oo.html) - O'Connor-Ott table chemical-potential conventions

## See Also

- Depends on: [Nucleon Blocking And EOS Conventions](nucleon-blocking-and-eos-conventions.md)
- Depends on: [StellarCollapse Table Adapter](../eos/stellarcollapse-table-adapter.md)
- See also: [Physics Model And Sources](physics-model-and-sources.md)
- See also: [CompOSE EOS Leakage Callback](compose-eos-leakage-callback.md)
- Validated by: [Neutrinos Tests And Reference Data](tests-and-reference-data.md)
