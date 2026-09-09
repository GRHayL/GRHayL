#include "ghl_m1.h"
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

static ghl_error_codes_t ghl_m1_neutrino_validate_lepton_weight(
      const ghl_m1_neutrino_rates *restrict rates) {

  switch(rates->species) {
    case ghl_m1_neutrino_nue:
      if(rates->lepton_weight != 1.0)
        return ghl_error_m1_microphysics_failure;
      return ghl_success;
    case ghl_m1_neutrino_anue:
      if(rates->lepton_weight != -1.0)
        return ghl_error_m1_microphysics_failure;
      return ghl_success;
    case ghl_m1_neutrino_nux:
      if(rates->lepton_weight != 0.0)
        return ghl_error_m1_microphysics_failure;
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

  if(rates == NULL || dYe_matter == NULL)
    return ghl_error_m1_null_pointer;

  if(!isfinite(dL_rad_cc) || !isfinite(baryon_density_conserved))
    return ghl_error_m1_invalid_state;

  if(baryon_density_conserved <= 0.0)
    return ghl_error_m1_invalid_state;

  ghl_error_codes_t error = ghl_m1_neutrino_validate_lepton_weight(rates);
  if(error != ghl_success)
    return error;

  if(!isfinite(rates->lepton_weight))
    return ghl_error_m1_microphysics_failure;

  const double candidate = -dL_rad_cc / baryon_density_conserved;
  if(!isfinite(candidate))
    return ghl_error_m1_invalid_state;

  *dYe_matter = candidate;
  return ghl_success;
}
