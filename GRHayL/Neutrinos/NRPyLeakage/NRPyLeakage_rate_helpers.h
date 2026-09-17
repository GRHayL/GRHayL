#ifndef NRPYLEAKAGE_RATE_HELPERS_H_
#define NRPYLEAKAGE_RATE_HELPERS_H_

#include <math.h>

#include "ghl_nrpyleakage.h"

/** Small positive fallback used by generated cgs leakage expressions. */
static const double nrpyl_cgs_numerical_floor = 1.0e-15;

/** Convert the cgs inverse-length floor to the public geometrized unit. */
static inline double nrpyl_opacity_floor_geom(void) {
  return NRPyLeakage_units_geom_to_cgs_L * nrpyl_cgs_numerical_floor;
}

/**
 * Compute the unsuppressed nucleon-bremsstrahlung number-emission rate.
 *
 * @param[in] T Temperature in MeV.
 * @param[in] rho_cgs Rest-mass density in g/cm^3.
 * @param[in] X_n Free-neutron mass fraction.
 * @param[in] X_p Free-proton mass fraction.
 * @return Number-emission rate in cgs units.
 */
static inline double nrpyl_bremsstrahlung_number_rate(
      const double T,
      const double rho_cgs,
      const double X_n,
      const double X_p) {
  const double composition
        = X_n * X_n + (28.0 / 3.0) * X_n * X_p + X_p * X_p;
  return NRPyLeakage_Brems_C1 * NRPyLeakage_Brems_zeta * pow(T, 4.5) * rho_cgs
         * rho_cgs * composition;
}

/**
 * Convert a bremsstrahlung number-emission rate to its energy moment.
 *
 * @param[in] T Temperature in MeV.
 * @param[in] number_rate Number-emission rate in cgs units.
 * @return Energy-emission rate in cgs units.
 */
static inline double nrpyl_bremsstrahlung_energy_rate(
      const double T,
      const double number_rate) {
  return NRPyLeakage_Brems_C2 * T * number_rate / NRPyLeakage_Brems_C1;
}

/**
 * Suppress a free emission rate by the leakage diffusion timescale.
 *
 * The cgs diffusion time is @f$t_{\rm diff}=6\tau^2/(c\kappa)@f$.
 * One free rate supplies both the numerator and its own inverse loss time;
 * this is especially important for the single-species heavy-lepton rate.
 *
 * @param[in] free_rate Optically thin emission rate.
 * @param[in] tau Optical depth for the corresponding moment and species.
 * @param[in] opacity Opacity for the corresponding moment and species.
 * @param[in] phase_space Equilibrium phase-space density.
 * @return Effective emission rate.
 */
static inline double nrpyl_effective_emission_rate(
      const double free_rate,
      const double tau,
      const double opacity,
      const double phase_space) {
  if(free_rate == 0.0)
    return 0.0;
  if(tau == 0.0)
    return free_rate;
  if(opacity == 0.0)
    return 0.0;

  const double diffusion_factor = 6.0 / NRPyLeakage_c_light;
  return free_rate
         / (tau * tau * diffusion_factor * free_rate
                  / (opacity * fmax(phase_space, 1.0000000000000001e-15))
            + 1.0);
}

/**
 * Combine the energy rates that cool the fluid.
 *
 * The heavy-lepton rate is per species, so the matter source includes all
 * four species: muon and tau neutrinos and their antiparticles.
 *
 * @param[in] Q_eff_nue Electron-neutrino effective energy rate.
 * @param[in] Q_eff_anue Electron-antineutrino effective energy rate.
 * @param[in] Q_eff_nux Single-species heavy-lepton effective energy rate.
 * @return Geometrized matter energy source.
 */
static inline double nrpyl_matter_energy_source(
      const double Q_eff_nue,
      const double Q_eff_anue,
      const double Q_eff_nux) {
  return NRPyLeakage_units_cgs_to_geom_Q
         * (-Q_eff_nue - Q_eff_anue - 4.0 * Q_eff_nux);
}

/** Replace a nonfinite opacity with the established positive floor. */
static inline double nrpyl_finite_opacity_or_floor(const double value) {
  return robust_isfinite(value) ? value : nrpyl_opacity_floor_geom();
}

/** Replace a nonfinite emitted or signed rate with the neutral value. */
static inline double nrpyl_finite_rate_or_zero(const double value) {
  return robust_isfinite(value) ? value : 0.0;
}

/** Apply the finite-output contract to every public opacity component. */
static inline bool nrpyl_sanitize_opacities(
  ghl_neutrino_opacities *restrict kappa) {
  bool replaced = false;
  for(int i = 0; i < 2; i++) {
    replaced |= !robust_isfinite(kappa->nue[i]);
    replaced |= !robust_isfinite(kappa->anue[i]);
    replaced |= !robust_isfinite(kappa->nux[i]);
    kappa->nue[i] = nrpyl_finite_opacity_or_floor(kappa->nue[i]);
    kappa->anue[i] = nrpyl_finite_opacity_or_floor(kappa->anue[i]);
    kappa->nux[i] = nrpyl_finite_opacity_or_floor(kappa->nux[i]);
  }
  return replaced;
}

/** Apply the finite-output contract to every public luminosity component. */
static inline bool nrpyl_sanitize_luminosities(
      ghl_neutrino_luminosities *restrict lum) {
  const bool replaced = !robust_isfinite(lum->nue)
                        || !robust_isfinite(lum->anue)
                        || !robust_isfinite(lum->nux);
  lum->nue = nrpyl_finite_rate_or_zero(lum->nue);
  lum->anue = nrpyl_finite_rate_or_zero(lum->anue);
  lum->nux = nrpyl_finite_rate_or_zero(lum->nux);
  return replaced;
}

/** Apply the finite-output contract to both public GRMHD source terms. */
static inline bool nrpyl_sanitize_sources(
      double *restrict R_source,
      double *restrict Q_source) {
  const bool replaced = !robust_isfinite(*R_source)
                        || !robust_isfinite(*Q_source);
  *R_source = nrpyl_finite_rate_or_zero(*R_source);
  *Q_source = nrpyl_finite_rate_or_zero(*Q_source);
  return replaced;
}

#endif // NRPYLEAKAGE_RATE_HELPERS_H_
