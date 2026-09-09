#include "ghl_m1_neutrino_implicit.h"

ghl_error_codes_t ghl_m1_neutrino_assemble_exchange(
      const ghl_m1_neutrino_state *restrict state_in,
      const ghl_m1_neutrino_state *restrict state_out,
      const ghl_m1_neutrino_rates *restrict rates,
      const double dL_rad_cc,
      const double sqrt_detgamma,
      const double baryon_density_conserved,
      ghl_m1_neutrino_exchange *restrict exchange) {

  if(state_in == NULL || state_out == NULL || rates == NULL || exchange == NULL)
    return ghl_error_m1_null_pointer;
  if(!isfinite(sqrt_detgamma) || sqrt_detgamma <= 0.0 ||
     !isfinite(dL_rad_cc))
    return ghl_error_m1_invalid_state;

  ghl_m1_neutrino_exchange candidate = {0};
  candidate.dN_rad_total = state_out->N - state_in->N;
  candidate.dL_rad_cc = dL_rad_cc;
  candidate.dE_rad = state_out->E - state_in->E;
  if(!isfinite(candidate.dN_rad_total) || !isfinite(candidate.dE_rad))
    return ghl_error_m1_invalid_state;
  for(int i = 0; i < 3; ++i) {
    candidate.dF_rad[i] = state_out->F[i] - state_in->F[i];
    if(!isfinite(candidate.dF_rad[i]))
      return ghl_error_m1_invalid_state;
  }

  candidate.dTau_matter = -sqrt_detgamma * candidate.dE_rad;
  if(!isfinite(candidate.dTau_matter))
    return ghl_error_m1_invalid_state;
  for(int i = 0; i < 3; ++i) {
    candidate.dS_matter[i] = -sqrt_detgamma * candidate.dF_rad[i];
    if(!isfinite(candidate.dS_matter[i]))
      return ghl_error_m1_invalid_state;
  }

  const ghl_error_codes_t error = ghl_m1_compute_neutrino_lepton_increment(
      rates, candidate.dL_rad_cc, baryon_density_conserved,
      &candidate.dYe_matter);
  if(error != ghl_success)
    return error;
  if(!isfinite(candidate.dYe_matter))
    return ghl_error_m1_invalid_state;

  *exchange = candidate;
  return ghl_success;
}
