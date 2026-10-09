#include "ghl_m1_neutrino_rate_backend.h"

#ifdef GHL_DISABLE_HDF5

#define BACKEND_STUB_UNUSED(value) ((void)(value))

ghl_error_codes_t ghl_m1_neutrino_rate_backend_initialize(void) {
  return ghl_error_used_disabled_hdf5;
}

ghl_error_codes_t ghl_m1_neutrino_rate_backend_temperature_from_eps(
      const ghl_eos_parameters *restrict eos,
      double rho,
      double Ye,
      double eps,
      double *restrict T) {
  BACKEND_STUB_UNUSED(eos);
  BACKEND_STUB_UNUSED(rho);
  BACKEND_STUB_UNUSED(Ye);
  BACKEND_STUB_UNUSED(eps);
  BACKEND_STUB_UNUSED(T);
  return ghl_error_used_disabled_hdf5;
}

ghl_error_codes_t ghl_m1_neutrino_rate_backend_thermo_from_T(
      const ghl_eos_parameters *restrict eos,
      double rho,
      double Ye,
      double T,
      double *restrict muhat,
      double *restrict mu_e,
      double *restrict mu_p,
      double *restrict mu_n,
      double *restrict X_n,
      double *restrict X_p) {
  BACKEND_STUB_UNUSED(eos);
  BACKEND_STUB_UNUSED(rho);
  BACKEND_STUB_UNUSED(Ye);
  BACKEND_STUB_UNUSED(T);
  BACKEND_STUB_UNUSED(muhat);
  BACKEND_STUB_UNUSED(mu_e);
  BACKEND_STUB_UNUSED(mu_p);
  BACKEND_STUB_UNUSED(mu_n);
  BACKEND_STUB_UNUSED(X_n);
  BACKEND_STUB_UNUSED(X_p);
  return ghl_error_used_disabled_hdf5;
}

ghl_error_codes_t ghl_m1_neutrino_rate_backend_assemble(
      int channel_mask,
      double nu_x_multiplicity,
      double rho,
      double T,
      double Ye,
      double muhat,
      double mu_e,
      double mu_p,
      double mu_n,
      double X_n,
      double X_p,
      ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count],
      double mismatch[ghl_m1_neutrino_species_count],
      bool mismatch_valid[ghl_m1_neutrino_species_count]) {
  BACKEND_STUB_UNUSED(channel_mask);
  BACKEND_STUB_UNUSED(nu_x_multiplicity);
  BACKEND_STUB_UNUSED(rho);
  BACKEND_STUB_UNUSED(T);
  BACKEND_STUB_UNUSED(Ye);
  BACKEND_STUB_UNUSED(muhat);
  BACKEND_STUB_UNUSED(mu_e);
  BACKEND_STUB_UNUSED(mu_p);
  BACKEND_STUB_UNUSED(mu_n);
  BACKEND_STUB_UNUSED(X_n);
  BACKEND_STUB_UNUSED(X_p);
  BACKEND_STUB_UNUSED(rates);
  BACKEND_STUB_UNUSED(mismatch);
  BACKEND_STUB_UNUSED(mismatch_valid);
  return ghl_error_used_disabled_hdf5;
}

#undef BACKEND_STUB_UNUSED

#else

#include "ghl_m1_nrpyleakage_kernel.h"
#include "ghl_nrpyeos_tabulated.h"
#include "ghl_nrpyleakage.h"
#include <float.h>
#include <math.h>
#include <string.h>

#define MEV_TO_ERG 1.602176634e-6

static double species_lepton_weight(const ghl_m1_neutrino_species_t species) {
  static const double weights[ghl_m1_neutrino_species_count]
        = { [ghl_m1_neutrino_nue] = 1.0,
            [ghl_m1_neutrino_anue] = -1.0,
            [ghl_m1_neutrino_nux] = 0.0 };
  return weights[species];
}

ghl_error_codes_t ghl_m1_neutrino_rate_backend_initialize(void) { return ghl_success; }

ghl_error_codes_t ghl_m1_neutrino_rate_backend_temperature_from_eps(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Ye,
      const double eps,
      double *restrict T) {
  double candidate = *T;
  const ghl_error_codes_t error
        = ghl_tabulated_compute_T_from_eps(eos, rho, Ye, eps, &candidate);
  if(error != ghl_success) {
    return error;
  }
  if(!isfinite(candidate) || candidate <= 0.0) {
    return ghl_error_m1_microphysics_failure;
  }
  *T = candidate;
  return ghl_success;
}

ghl_error_codes_t ghl_m1_neutrino_rate_backend_thermo_from_T(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Ye,
      const double T,
      double *restrict muhat,
      double *restrict mu_e,
      double *restrict mu_p,
      double *restrict mu_n,
      double *restrict X_n,
      double *restrict X_p) {
  double candidate_muhat = 0.0, candidate_mu_e = 0.0, candidate_mu_p = 0.0;
  double candidate_mu_n = 0.0, candidate_X_n = 0.0, candidate_X_p = 0.0;
  ghl_error_codes_t error = ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T(
        eos, rho, Ye, T, &candidate_muhat, &candidate_mu_e, &candidate_mu_p,
        &candidate_mu_n, &candidate_X_n, &candidate_X_p);
  if(error != ghl_success) {
    return error;
  }
  if(!isfinite(candidate_mu_e) || !isfinite(candidate_mu_p) || !isfinite(candidate_mu_n)
     || !isfinite(candidate_muhat)) {
    return ghl_error_m1_microphysics_failure;
  }
  double normalized_X_n = candidate_X_n, normalized_X_p = candidate_X_p;
  error = ghl_nrpyleakage_normalize_nucleon_fractions(
        candidate_X_n, candidate_X_p, &normalized_X_n, &normalized_X_p);
  if(error != ghl_success) {
    return ghl_error_m1_microphysics_failure;
  }
  *muhat = candidate_muhat;
  *mu_e = candidate_mu_e;
  *mu_p = candidate_mu_p;
  *mu_n = candidate_mu_n;
  *X_n = normalized_X_n;
  *X_p = normalized_X_p;
  return ghl_success;
}

ghl_error_codes_t ghl_m1_neutrino_rate_backend_assemble(
      const int channel_mask,
      const double nu_x_multiplicity,
      const double rho,
      const double T,
      const double Ye,
      const double muhat,
      const double mu_e,
      const double mu_p,
      const double mu_n,
      const double X_n,
      const double X_p,
      ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count],
      double beta_kirchhoff_relative_mismatch[ghl_m1_neutrino_species_count],
      bool beta_kirchhoff_mismatch_valid[ghl_m1_neutrino_species_count]) {

  ghl_m1_nrpyleakage_thermo_state thermo;
  ghl_error_codes_t err = ghl_m1_nrpyleakage_build_thermo_state_from_eos_quantities(
        rho, Ye, T, muhat, mu_e, mu_p, mu_n, X_n, X_p, &thermo);
  if(err != ghl_success) {
    return err;
  }

  const double eta[ghl_m1_nrpyleakage_species_count]
        = { (mu_e - muhat) / T, -(mu_e - muhat) / T, 0.0 };
  ghl_m1_nrpyleakage_raw_rates raw;
  err = ghl_m1_nrpyleakage_compute_raw_rates_from_thermo_with_mask(
        &thermo, eta, channel_mask, &raw);
  if(err != ghl_success) {
    return err;
  }
  /* The raw bridge exposes one unsummed heavy species.  The provider owns
   * the sole conversion to the configured four-flavor aggregate below. */

  const double L0 = NRPyLeakage_units_geom_to_cgs_L;
  const double t0 = NRPyLeakage_units_geom_to_cgs_T;
  const double E0
        = NRPyLeakage_units_geom_to_cgs_M * NRPyLeakage_c_light * NRPyLeakage_c_light;
  const double volume = L0 * L0 * L0;
  const double number_emissivity_conversion = volume * t0;
  const double energy_density_conversion = MEV_TO_ERG * volume / E0;
  const double energy_emissivity_conversion = MEV_TO_ERG * volume * t0 / E0;
  const int mask = channel_mask;

  for(int s = 0; s < ghl_m1_neutrino_species_count; ++s) {
    const ghl_m1_nrpyleakage_species_raw_rates *const rr = &raw.species[s];
    ghl_m1_neutrino_rates *const out = &rates[s];
    const double multiplicity = s == ghl_m1_neutrino_nux ? nu_x_multiplicity : 1.0;
    beta_kirchhoff_relative_mismatch[s] = 0.0;
    beta_kirchhoff_mismatch_valid[s] = false;
    *out = (ghl_m1_neutrino_rates){ 0 };
    out->species = (ghl_m1_neutrino_species_t)s;
    out->lepton_weight = species_lepton_weight(out->species);
    const double raw_n_eq = multiplicity * rr->n_eq_cgs * volume;
    const double raw_J_eq = multiplicity * rr->J_eq_mev_cgs * energy_density_conversion;
    /* Successful raw evaluation implies T < 1e62 from the T^5 product and
     * |eta| < 746 from both electron-species exponential tails. The FD
     * approximations then give F2 < 1.4e8 and F3 < 8e10, bounding raw n by
     * 1e226 and raw J by 1e291. The fixed four-flavor multiplicity and unit
     * conversions cannot overflow either target; underflow remains possible. */
    if(raw_n_eq <= 0.0 || raw_J_eq <= 0.0) {
      return ghl_error_m1_microphysics_failure;
    }
    /* Preserve the two raw FD moments.  The mean is derived from the same
     * converted targets so validation and detailed balance use one physical
     * equilibrium state.  A zero/overflowing target is a genuine binary64
     * representability failure, not a reason to substitute an energy floor. */
    out->n_eq = raw_n_eq;
    out->J_eq = raw_J_eq;
    out->mean_energy = out->J_eq / out->n_eq;
    /* The same FD bounds make this mean finite. Positive converted J also
     * requires nonzero T^4, hence T > 1e-81; with F3/F2 of order one or
     * larger and the fixed energy conversion, the ratio cannot underflow. */

    if(s != ghl_m1_neutrino_nux) {
      const double beta = rr->eta_N_beta_cgs * number_emissivity_conversion;
      const double kirchhoff = rr->kappa_a_N_cc_cgs * L0 * out->n_eq;
      beta_kirchhoff_relative_mismatch[s]
            = fabs(beta - kirchhoff) / fmax(fmax(fabs(beta), fabs(kirchhoff)), DBL_MIN);
      beta_kirchhoff_mismatch_valid[s] = true;
    }

    double kappa_a_N = 0.0;
    double kappa_a_E = 0.0;
    if((mask & ghl_neutrino_rate_channel_charged_current) != 0
       && s != ghl_m1_neutrino_nux) {
      out->kappa_a_N_cc = rr->kappa_a_N_cc_cgs * L0;
      kappa_a_N += out->kappa_a_N_cc;
      kappa_a_E += rr->kappa_a_E_cc_cgs * L0;
    }

    if((mask & ghl_neutrino_rate_channel_pair) != 0 && s == ghl_m1_neutrino_nux) {
      const double eta_N_channel
            = multiplicity * rr->eta_N_pair_cgs * number_emissivity_conversion;
      const double eta_E_channel
            = multiplicity * rr->eta_E_pair_mev_cgs * energy_emissivity_conversion;
      kappa_a_N += eta_N_channel / out->n_eq;
      kappa_a_E += eta_E_channel / out->J_eq;
    }
    if((mask & ghl_neutrino_rate_channel_plasmon) != 0 && s == ghl_m1_neutrino_nux) {
      const double eta_N_channel
            = multiplicity * rr->eta_N_plasmon_cgs * number_emissivity_conversion;
      const double eta_E_channel
            = multiplicity * rr->eta_E_plasmon_mev_cgs * energy_emissivity_conversion;
      kappa_a_N += eta_N_channel / out->n_eq;
      kappa_a_E += eta_E_channel / out->J_eq;
    }
    if((mask & ghl_neutrino_rate_channel_bremsstrahlung) != 0
       && s == ghl_m1_neutrino_nux) {
      const double eta_N_channel
            = multiplicity * rr->eta_N_brems_cgs * number_emissivity_conversion;
      const double eta_E_channel
            = multiplicity * rr->eta_E_brems_mev_cgs * energy_emissivity_conversion;
      kappa_a_N += eta_N_channel / out->n_eq;
      kappa_a_E += eta_E_channel / out->J_eq;
    }

    /* Electron-flavor emission remains separate from the aggregate
     * charged-current/scattering coefficients.  The raw adapter supplies a
     * common number emissivity for the electron pair; copy both raw number
     * and energy rates without rebuilding either from a grey opacity. */
    if(s != ghl_m1_neutrino_nux) {
      if((mask & ghl_neutrino_rate_channel_pair) != 0) {
        out->eta_N_pair[ghl_m1_neutrino_pair_process_pair]
              = rr->eta_N_pair_cgs * number_emissivity_conversion;
        out->eta_E_pair[ghl_m1_neutrino_pair_process_pair]
              = rr->eta_E_pair_mev_cgs * energy_emissivity_conversion;
      }
      if((mask & ghl_neutrino_rate_channel_plasmon) != 0) {
        out->eta_N_pair[ghl_m1_neutrino_pair_process_plasmon]
              = rr->eta_N_plasmon_cgs * number_emissivity_conversion;
        out->eta_E_pair[ghl_m1_neutrino_pair_process_plasmon]
              = rr->eta_E_plasmon_mev_cgs * energy_emissivity_conversion;
      }
      if((mask & ghl_neutrino_rate_channel_bremsstrahlung) != 0) {
        out->eta_N_pair[ghl_m1_neutrino_pair_process_bremsstrahlung]
              = rr->eta_N_brems_cgs * number_emissivity_conversion;
        out->eta_E_pair[ghl_m1_neutrino_pair_process_bremsstrahlung]
              = rr->eta_E_brems_mev_cgs * energy_emissivity_conversion;
      }
    }

    out->kappa_a_N = kappa_a_N;
    out->kappa_a_E = kappa_a_E;
    if((mask & ghl_neutrino_rate_channel_nucleon_scattering) != 0) {
      out->kappa_s = L0 * (rr->kappa_s_E_neutron_cgs + rr->kappa_s_E_proton_cgs);
    }
    out->kappa_tr = out->kappa_a_E + out->kappa_s;
    out->eta_N = out->kappa_a_N * out->n_eq;
    out->eta_E = out->kappa_a_E * out->J_eq;
    out->eta_N_cc = out->kappa_a_N_cc * out->n_eq;
  }
  return ghl_success;
}

#endif /* GHL_DISABLE_HDF5 */
