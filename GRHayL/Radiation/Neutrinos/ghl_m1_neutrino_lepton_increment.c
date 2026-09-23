#include "ghl_m1.h"
#include "ghl_m1_neutrino_implicit.h"
#include <float.h>

/*
 * Matter Delta Y_e recommendation. dL_rad_cc is already
 * species-signed:
 *   Delta Y_e_matter = -Delta L_rad_cc / baryon_density_conserved
 *
 * Sign consequences:
 *   - emitting nu_e     (dL_rad_cc > 0, lw = +1) -> Delta Y_e < 0
 *   - absorbing nu_e    (dL_rad_cc < 0, lw = +1) -> Delta Y_e > 0
 *   - emitting nu_bar_e (dL_rad_cc < 0) -> Delta Y_e > 0
 *   - absorbing nu_bar_e(dL_rad_cc > 0) -> Delta Y_e < 0
 *   - nu_x (lw = 0) gives Delta Y_e = 0 when dL_rad_cc = 0; the caller
 *     must supply zero charged-current exchange for this species.
 *
 * The provider supplies the lepton weight via ghl_m1_neutrino_rates and is
 * responsible for ensuring species consistency. baryon_density_conserved
 * must be positive; the host may substitute its own densitized baryon-mass
 * normalization when applying the increment.
 *
 */

static ghl_error_codes_t
ghl_m1_neutrino_validate_lepton_weight(const ghl_m1_neutrino_rates *restrict rates) {

  switch(rates->species) {
    case ghl_m1_neutrino_nue:
      if(rates->lepton_weight != 1.0) {
        return ghl_error_m1_microphysics_failure;
      }
      return ghl_success;
    case ghl_m1_neutrino_anue:
      if(rates->lepton_weight != -1.0) {
        return ghl_error_m1_microphysics_failure;
      }
      return ghl_success;
    case ghl_m1_neutrino_nux:
      if(rates->lepton_weight != 0.0) {
        return ghl_error_m1_microphysics_failure;
      }
      return ghl_success;
    default:
      return ghl_error_m1_microphysics_failure;
  }
}

ghl_error_codes_t ghl_m1_compute_neutrino_lepton_increment(
      const ghl_m1_neutrino_rates *restrict rates,
      const double dL_rad_cc,
      const double baryon_density_conserved,
      double *restrict dYe_matter) {

  if(rates == NULL || dYe_matter == NULL) {
    return ghl_error_m1_null_pointer;
  }

  if(!isfinite(dL_rad_cc) || !isfinite(baryon_density_conserved)) {
    return ghl_error_m1_invalid_state;
  }

  if(baryon_density_conserved <= 0.0) {
    return ghl_error_m1_invalid_state;
  }

  ghl_error_codes_t error = ghl_m1_neutrino_validate_lepton_weight(rates);
  if(error != ghl_success) {
    return error;
  }

  if(!isfinite(rates->lepton_weight)) {
    return ghl_error_m1_microphysics_failure;
  }

  /* Heavy flavor carries no electron-lepton number.  The documented contract
   * above requires the caller to supply zero charged-current exchange for this
   * species; enforce it here so a malformed bundle cannot silently publish a
   * nonzero composition source.  Signed zero satisfies the check. */
  if(rates->species == ghl_m1_neutrino_nux && dL_rad_cc != 0.0) {
    return ghl_error_m1_microphysics_failure;
  }

  const double candidate = -dL_rad_cc / baryon_density_conserved;
  if(!isfinite(candidate)) {
    return ghl_error_m1_invalid_state;
  }

  *dYe_matter = candidate;
  return ghl_success;
}

ghl_error_codes_t ghl_m1_neutrino_charged_current_lepton_delta(
      const ghl_m1_neutrino_rates *restrict rates,
      const double dt_alpha,
      const double number_initial,
      const bool number_projected,
      const double physical_number_endpoint,
      const double physical_number_gamma,
      double *restrict dL_rad_cc) {

  if(rates == NULL || dL_rad_cc == NULL) {
    return ghl_error_m1_null_pointer;
  }
  if(!isfinite(number_initial) || number_initial < 0.0
     || !isfinite(dt_alpha) || dt_alpha < 0.0 || !isfinite(physical_number_endpoint)
     || physical_number_endpoint < 0.0 || !isfinite(physical_number_gamma)
     || physical_number_gamma <= 0.0) {
    return ghl_error_m1_invalid_state;
  }

  /* Validated single-species electron rates contain only charged-current
   * number reactions. Ordinary backward Euler therefore transfers exactly
   * the un-repaired change in N. Re-evaluating emission minus absorption at
   * its rounded endpoint loses that change in the stiff equilibrium limit.
   * Heavy flavor has no charged-current exchange. The optional mean-energy
   * projection is not backward Euler and retains its endpoint-source rule. */
  double dN_cc = 0.0;
  if(rates->lepton_weight != 0.0) {
    if(number_projected) {
      dN_cc = dt_alpha
              * (rates->eta_N_cc
                 - rates->kappa_a_N_cc * physical_number_endpoint / physical_number_gamma);
    }
    else {
      dN_cc = physical_number_endpoint - number_initial;
    }
  }
  const double candidate = rates->lepton_weight * dN_cc;
  if(!isfinite(candidate)) {
    return ghl_error_m1_invalid_state;
  }

  *dL_rad_cc = candidate;
  return ghl_success;
}
