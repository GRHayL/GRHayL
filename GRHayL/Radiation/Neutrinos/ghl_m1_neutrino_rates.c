#include "ghl_m1.h"
#include "ghl_m1_neutrino_implicit.h"
#include "../ghl_m1_utils.h"
#include <float.h>
#include <fenv.h>

/*
 * Validation of a provider-supplied frozen
 * ghl_m1_neutrino_rates bundle.
 *
 * Validation rules:
 *   - all fields finite;
 *   - emissivities, opacities, equilibrium values >= 0; mean_energy > 0;
 *   - species in {nue, anue, nux};
 *   - lepton_weight matches the species exactly:
 *       nu_e   -> +1
 *       nu_bar_e -> -1
 *       nu_x   ->  0
 *   - total and charged-current-subset Kirchhoff identities, transport-opacity
 *     sum, and J_eq = n_eq*mean_energy use the specified binary64 comparison;
 *
 * On validation failure with non-NULL diagnostics, the provider_validation_failures
 * counter is incremented and ghl_error_m1_microphysics_failure is returned.
 * The bundle is never modified.
 */

static ghl_error_codes_t ghl_m1_neutrino_check_lepton_weight(
      const ghl_m1_neutrino_species_t species,
      const double lepton_weight) {

  switch(species) {
    case ghl_m1_neutrino_nue:
      if(lepton_weight != 1.0)
        return ghl_error_m1_microphysics_failure;
      return ghl_success;
    case ghl_m1_neutrino_anue:
      if(lepton_weight != -1.0)
        return ghl_error_m1_microphysics_failure;
      return ghl_success;
    case ghl_m1_neutrino_nux:
      if(lepton_weight != 0.0)
        return ghl_error_m1_microphysics_failure;
      return ghl_success;
    default:
      return ghl_error_m1_microphysics_failure;
  }
}

static bool rate_close(const double x, const double y) {
  return fabs(x-y) <= 128.0*DBL_EPSILON*ghl_m1_max(fabs(x), fabs(y));
}

static bool rate_expected_product(
      const double supplied, const double a, const double b,
      bool *restrict representational_underflow) {
  *representational_underflow = false;
  if(a == 0.0 || b == 0.0)
    return supplied == 0.0;
  volatile double expected_storage = a*b;
  const double expected = expected_storage;
  if(!isfinite(expected))
    return false;
  if(expected == 0.0) {
    *representational_underflow = true;
    return supplied == 0.0;
  }
  return rate_close(supplied, expected);
}

static bool rate_expected_sum(
      const double supplied, const double a, const double b) {
  volatile double expected_storage = a+b;
  const double expected = expected_storage;
  if(!isfinite(expected))
    return false;
  if(a == 0.0 && b == 0.0)
    return supplied == 0.0;
  return rate_close(supplied, expected);
}

ghl_error_codes_t ghl_m1_validate_neutrino_rates(
      const ghl_m1_neutrino_rates *restrict rates,
      ghl_m1_neutrino_diagnostics *restrict diagnostics) {

  if(rates == NULL)
    return ghl_error_m1_null_pointer;

  /* Species range */
  if(rates->species != ghl_m1_neutrino_nue &&
     rates->species != ghl_m1_neutrino_anue &&
     rates->species != ghl_m1_neutrino_nux) {
    if(diagnostics != NULL)
      diagnostics->provider_validation_failures++;
    return ghl_error_m1_microphysics_failure;
  }

  /* Finiteness of all rate fields */
  if(!isfinite(rates->eta_N)        || !isfinite(rates->eta_E)        ||
     !isfinite(rates->kappa_a_N)    || !isfinite(rates->kappa_a_E)    ||
     !isfinite(rates->kappa_s)      || !isfinite(rates->kappa_tr)     ||
     !isfinite(rates->n_eq)         || !isfinite(rates->J_eq)         ||
     !isfinite(rates->mean_energy)  || !isfinite(rates->lepton_weight) ||
     !isfinite(rates->eta_N_cc) || !isfinite(rates->kappa_a_N_cc)) {
    if(diagnostics != NULL)
      diagnostics->provider_validation_failures++;
    return ghl_error_m1_microphysics_failure;
  }

  /* Nonnegativity of emissivities and equilibrium targets */
  if(rates->eta_N  < 0.0 || rates->eta_E  < 0.0 ||
     rates->n_eq   < 0.0 || rates->J_eq   < 0.0) {
    if(diagnostics != NULL)
      diagnostics->provider_validation_failures++;
    return ghl_error_m1_microphysics_failure;
  }

  /* Nonnegativity of opacities */
  if(rates->kappa_a_N < 0.0 || rates->kappa_a_E < 0.0 ||
     rates->kappa_s   < 0.0 || rates->kappa_tr  < 0.0 ||
     rates->eta_N_cc < 0.0 || rates->kappa_a_N_cc < 0.0) {
    if(diagnostics != NULL)
      diagnostics->provider_validation_failures++;
    return ghl_error_m1_microphysics_failure;
  }

  /* Mean energy strictly positive */
  if(rates->mean_energy <= 0.0) {
    if(diagnostics != NULL)
      diagnostics->provider_validation_failures++;
    return ghl_error_m1_microphysics_failure;
  }

  for(int process = 0; process < ghl_m1_neutrino_pair_process_count; ++process) {
    const double eta_N = rates->eta_N_pair[process];
    const double eta_E = rates->eta_E_pair[process];
    if(!isfinite(eta_N) || !isfinite(eta_E) || eta_N < 0.0 || eta_E < 0.0 ||
       ((eta_N != 0.0 || eta_E != 0.0) &&
        (rates->species == ghl_m1_neutrino_nux ||
         rates->n_eq <= 0.0 || rates->J_eq <= 0.0))) {
      if(diagnostics != NULL) {
        diagnostics->provider_validation_failures++;
      }
      return ghl_error_m1_microphysics_failure;
    }
  }

  /* Species / lepton-weight consistency */
  if(ghl_m1_neutrino_check_lepton_weight(rates->species,
                                        rates->lepton_weight)
     != ghl_success) {
    if(diagnostics != NULL)
      diagnostics->provider_validation_failures++;
    return ghl_error_m1_microphysics_failure;
  }

  bool product_underflow[4] = {false, false, false, false};
  if(fegetround() != FE_TONEAREST ||
     rates->eta_N_cc > rates->eta_N ||
     rates->kappa_a_N_cc > rates->kappa_a_N ||
     (rates->species == ghl_m1_neutrino_nux &&
      (rates->eta_N_cc != 0.0 || rates->kappa_a_N_cc != 0.0)) ||
     !rate_expected_sum(rates->kappa_tr, rates->kappa_a_E, rates->kappa_s) ||
     !rate_expected_product(rates->eta_N, rates->kappa_a_N, rates->n_eq,
                            &product_underflow[0]) ||
     !rate_expected_product(rates->eta_E, rates->kappa_a_E, rates->J_eq,
                            &product_underflow[1]) ||
     !rate_expected_product(rates->eta_N_cc, rates->kappa_a_N_cc, rates->n_eq,
                            &product_underflow[2]) ||
     !rate_expected_product(rates->J_eq, rates->n_eq, rates->mean_energy,
                            &product_underflow[3])) {
    if(diagnostics != NULL)
      diagnostics->provider_validation_failures++;
    return ghl_error_m1_microphysics_failure;
  }

  if(diagnostics != NULL) {
    for(int i = 0; i < 4; i++)
      diagnostics->rate_product_underflows += product_underflow[i];
  }

  return ghl_success;
}

ghl_error_codes_t ghl_m1_neutrino_validate_single_species_rates(
      const ghl_m1_neutrino_rates *restrict rates,
      ghl_m1_neutrino_diagnostics *restrict diagnostics) {

  const ghl_error_codes_t error = ghl_m1_validate_neutrino_rates(rates, diagnostics);
  if(error != ghl_success) {
    return error;
  }
  if(rates->species == ghl_m1_neutrino_nux) {
    return ghl_success;
  }

  /* Neither separated pair coefficients nor legacy non-CC electron number
   * rates can be evaluated from one species alone. */
  bool needs_partner = rates->eta_N != rates->eta_N_cc ||
                       rates->kappa_a_N != rates->kappa_a_N_cc;
  for(int process = 0; process < ghl_m1_neutrino_pair_process_count; ++process) {
    needs_partner |= rates->eta_N_pair[process] != 0.0 ||
                     rates->eta_E_pair[process] != 0.0;
  }
  if(needs_partner) {
    if(diagnostics != NULL) {
      diagnostics->provider_validation_failures++;
    }
    return ghl_error_m1_microphysics_failure;
  }
  return ghl_success;
}
