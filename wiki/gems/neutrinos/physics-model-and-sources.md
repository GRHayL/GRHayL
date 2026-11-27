# Neutrinos Physics Model And Sources

> Grey neutrino-leakage reactions, species and `[0]`/`[1]` slot meanings, effective rates, Valencia source projections, and cgs suppression-ratio units. · Status: confirmed
> Up: [Neutrinos](index.md)

## Summary

NRPyLeakage implements a grey neutrino-leakage model that computes local free emission rates, evaluates energy-averaged transport opacities and optical depths, and suppresses free rates by local diffusion-to-loss-time ratios. It is not a Boltzmann, Monte Carlo, or moment-transport solver and includes neither nonlocal absorption heating nor neutrino momentum deposition in optically thin ejecta. Current GRHayL C source, headers, and tests are the implementation authority; the ancestral NRPy+ notebooks and accompanying papers provide derivation background.

## Detail

### Physical Model And Species Convention

Local free emission includes electron and positron capture on nucleons (beta processes), electron-positron pair annihilation, transverse plasmon decay, and nucleon-nucleon bremsstrahlung. Transport opacities include neutrino absorption on nucleons for $\nu_e$ and $\bar\nu_e$, plus neutrino scattering on free neutrons and protons for all species; heavy-lepton species have scattering opacity without beta absorption.

Public containers order species as `nue` ($\nu_e$), `anue` ($\bar\nu_e$), and `nux` (one representative heavy-lepton species out of $\nu_\mu, \bar\nu_\mu, \nu_\tau, \bar\nu_\tau$). The formulation assumes all four heavy-lepton species contribute equally: `nux` rates and `lum.nux` represent one species, `Q_source` multiplies the effective `nux` cooling rate by four, and heavy-lepton number rates do not enter `R_source`.

Pointwise luminosities returned in `ghl_neutrino_luminosities` equal `W * alpha^2 * sqrt(det(gamma)) * Q_eff` converted to geometric units. The routine takes no cell volume and performs no spatial sum; callers integrate over coordinate volume without multiplying by a second `sqrt(gamma)`.

Evidence (descriptive behavior): `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_luminosities.c::NRPyLeakage_compute_neutrino_luminosities`; corroboration: `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c`.

### Number And Energy Slots And Optical Depths

Although `ghl_neutrino_opacities` leaves `[0]` and `[1]` unnamed, the implementation assigns:

- `[0]` ($j=0$, neutrino-number transport): used for `nue` and `anue` equilibrium-degeneracy interpolation, number diffusion time, and `R_source`. `nux[0]` opacity is written and its optical depth updated, though current source and luminosity routines do not consume `nux[0]`.
- `[1]` ($j=1$, neutrino-energy transport): used for energy diffusion time, `Q_source`, and luminosity suppression.

Slots `[0]` and `[1]` use distinct Fermi-Dirac energy moments and cannot be interchanged. Schematically, total transport opacity (inverse mean free path) is:

$$\kappa_{t,j}^{\nu_e} = \kappa_{s,j}^{\nu_e n} + \kappa_{s,j}^{\nu_e p} + \kappa_{a,j}^{\nu_e n}, \qquad \kappa_{t,j}^{\bar\nu_e} = \kappa_{s,j}^{\bar\nu_e n} + \kappa_{s,j}^{\bar\nu_e p} + \kappa_{a,j}^{\bar\nu_e p}, \qquad \kappa_{t,j}^{\nu_x} = \kappa_{s,j}^{\nu_x n} + \kappa_{s,j}^{\nu_x p}.$$

Dimensionless optical depth represents $\tau_j^{\nu_i} = \int ds\,\kappa_{t,j}^{\nu_i}$. `NRPyLeakage_optical_depths_PathOfLeastResistance` performs one six-face-neighbor update with no diagonal paths; global iteration to convergence, outer-boundary initialization, and AMR synchronization remain caller responsibilities.

Evidence (descriptive behavior): `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_optical_depths_PathOfLeastResistance.c::NRPyLeakage_optical_depths_PathOfLeastResistance`; corroboration: `Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c::run_unit_test`.

### Effective Rates And Heavy-Lepton Timescale

For every species and transport moment, leakage suppresses the free rate by the local diffusion-to-loss-time ratio:

$$R_{\nu_i}^{\rm eff} = \frac{R_{\nu_i}^{\rm free}}{1 + t_{\nu_i,0}^{\rm diff} / t_{\nu_i,0}^{\rm loss}}, \qquad Q_{\nu_i}^{\rm eff} = \frac{Q_{\nu_i}^{\rm free}}{1 + t_{\nu_i,1}^{\rm diff} / t_{\nu_i,1}^{\rm loss}}.$$

In cgs, $t_{\nu_i,j}^{\rm diff} = 6(\tau_j^{\nu_i})^2 / (c\,\kappa_{t,j}^{\nu_i})$, while inverse loss times are $1 / t_{\nu_i,0}^{\rm loss} = R_{\nu_i}^{\rm free} / n_{\nu_i}$ and $1 / t_{\nu_i,1}^{\rm loss} = Q_{\nu_i}^{\rm free} / e_{\nu_i}$ using equilibrium neutrino number and energy densities $n_{\nu_i}$ and $e_{\nu_i}$. Small $\tau$ recovers free streaming; large $\tau$ gives diffusion-limited leakage.

Both the matter-source and luminosity routines evaluate the single-species heavy-lepton free energy rate `Q_free_nux` once and use it both as the `nux` effective-energy numerator and in its inverse loss time $(t_{\nu_x,1}^{\rm loss})^{-1} = Q_{\nu_x}^{\rm free} / e_{\nu_x}$ paired with `nux` energy optical depth and opacity. Luminosity reference generation exercises nonzero optical depths in `[1, 1000]` across this suppression path, whereas the optically thin gas test sets all optical depths to zero. Electron-neutrino and antineutrino degeneracies interpolate between zero (transparent) and equilibrium using `exp(-tau.nue[0])` and `exp(-tau.anue[0])`; heavy-lepton degeneracy is zero.

Evidence (descriptive behavior): `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_rate_helpers.h::nrpyl_effective_emission_rate`; corroboration: `Unit_Tests/unit_test_nrpyleakage_luminosities.c::generate_test_data`.

### Matter-Source Signs And Valencia Projections

Before baryon-mass and geometric-unit scaling, the fluid-frame source convention is:

$$R_{\rm source} = R_{\bar\nu_e}^{\rm eff} - R_{\nu_e}^{\rm eff}, \qquad Q_{\rm source} = -\left(Q_{\nu_e}^{\rm eff} + Q_{\bar\nu_e}^{\rm eff} + 4Q_{\nu_x}^{\rm eff}\right).$$

Positive `R_source` raises matter $Y_e$; positive cooling yields negative `Q_source`, reducing internal energy. The C implementation multiplies `R_source` by the atomic mass unit and cgs-to-geometric conversion and converts `Q_source` to geometric units. Under $\rho_b = m_u n_b$, `R_source` sources $\rho_b Y_e$ directly (an EOS using a different baryon-mass reference requires consistent adaptation). In the homogeneous optically thin test, $dY_e/dt = R_{\rm source}/\rho$ and $d\epsilon/dt = Q_{\rm source}/\rho$.

For the isotropic fluid-frame radiation closure $G^\mu = Q_{\rm source} u^\mu$ (Siegel & Metzger, Eqs. (1)-(3), (7), (19)), additions to the densitized Valencia conservative RHS are:

$$\Delta_D = 0, \qquad \Delta_\tau = \alpha\sqrt{\gamma}\,W Q_{\rm source}, \qquad \Delta_{S_j} = \alpha\sqrt{\gamma}\,u_j Q_{\rm source}, \qquad \Delta_{DY_e} = \alpha\sqrt{\gamma}\,R_{\rm source},$$

where $D = \rho_b W$, $W = \alpha u^0$, and $u_j = W v_j$ for the covariant Eulerian Valencia velocity. Because GRHayL stores transport velocity `prims.vU[i] = u^i / u^0`, form the contravariant Eulerian velocity as `(prims.vU[i] + beta^i) / alpha` before lowering its index. Returned scalars carry no metric factors, velocity factors, or timestep, and already include baryon mass and the factor of four for heavy leptons.

Evidence (scientific/interface requirement): `GRHayL/include/ghl_nrpyleakage.h::NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms`; corroboration: `Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c::run_unit_test`.

### Units, Suppression-Ratio Scaling, And Model Limits

Leakage constants assume geometrized units ($G = c = M_\odot = 1$) in `GRHayL/include/ghl_nrpyleakage.h`. These macros are not numerically identical to the `CODE_TO_CGS_*` / `CGS_TO_CODE_*` macros in `GRHayL/include/ghl_nrpyeos_tabulated.h`; ports must update density, length, number-rate, and energy-rate conversions as one coherent set rather than mixing EOS and leakage macros.

The diffusion-to-loss-time ratio is built entirely in cgs using `6.0 / NRPyLeakage_c_light`, so $t^{\rm diff} = 6\tau^2 / (c\,\kappa_{\rm cgs})$ multiplies an inverse loss rate in $\mathrm{s}^{-1}$ to form a dimensionless ratio; opacity is converted to inverse geometric length only when writing `kappa`. Applying `NRPyLeakage_units_geom_to_cgs_T` or `NRPyLeakage_units_geom_to_cgs_L` inconsistently to this prefactor introduces an error of `NRPyLeakage_units_cgs_to_geom_T` ($\approx 2.03 \times 10^5$).

Model limits: grey energy-averaged rates only; no nonlocal heating or momentum deposition; four equal heavy-lepton species; limited reaction and opacity channels; focused blocking input checks rather than full thermodynamic validation; `EnsureFinite` guards are numerical fallbacks rather than EOS-convention checks; and six-face-neighbor path updates without built-in global convergence or outer boundaries.

Evidence (descriptive behavior): `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_rate_helpers.h::NRPyLeakage_c_light`; corroboration: `GRHayL/include/ghl_nrpyleakage.h::NRPyLeakage_units_cgs_to_geom_T`.

## Sources

- [GRHayL/include/ghl_radiation.h](../../../GRHayL/include/ghl_radiation.h) - `ghl_neutrino_luminosities`, `ghl_neutrino_opacities`, `ghl_neutrino_optical_depths`
- [GRHayL/include/ghl_nrpyleakage.h](../../../GRHayL/include/ghl_nrpyleakage.h) - `NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms`, `NRPyLeakage_c_light`, `NRPyLeakage_units_cgs_to_geom_T`
- [GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_rate_helpers.h](../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_rate_helpers.h) - `nrpyl_effective_emission_rate`, cgs diffusion-time prefactor
- [GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c](../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c) - `R_source`, `Q_source`, `[0]`/`[1]` opacity writeback
- [GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_luminosities.c](../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_luminosities.c) - pointwise luminosity scaling and suppression
- [GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_optical_depths_PathOfLeastResistance.c](../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_optical_depths_PathOfLeastResistance.c) - six-face minimum optical-depth update
- [Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c](../../../Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c) - `run_unit_test`
- [Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c](../../../Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c) - `run_unit_test`
- [Unit_Tests/unit_test_nrpyleakage_luminosities.c](../../../Unit_Tests/unit_test_nrpyleakage_luminosities.c) - `generate_test_data`
- [Siegel & Metzger (arXiv:1711.00868)](https://arxiv.org/abs/1711.00868) - Eqs. (1)-(3), (7), (19) Valencia leakage closure
- [Werneck et al. (arXiv:2208.14487)](https://arxiv.org/abs/2208.14487) - NRPyLeakage formulation
- [Werneck et al. PRD 107, 044037](https://doi.org/10.1103/PhysRevD.107.044037) - published NRPyLeakage reference
- [Ruffert, Janka & Schaefer (1996)](https://ui.adsabs.harvard.edu/abs/1996A%26A...311..532R/abstract) - grey leakage background

## See Also

- Depends on: [EOS Inputs And Chemical Potentials](eos-inputs-and-chemical-potentials.md)
- Depends on: [Nucleon Blocking And EOS Conventions](nucleon-blocking-and-eos-conventions.md)
- Implements: [Implementation Flow](implementation-flow.md)
- Validated by: [Neutrinos Tests And Reference Data](tests-and-reference-data.md)
- See also: [Variables And Conventions](../../physics/variables-and-conventions.md)
