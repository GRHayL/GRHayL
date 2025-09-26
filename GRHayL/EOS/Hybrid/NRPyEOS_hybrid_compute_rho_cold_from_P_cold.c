#include "ghl_nrpyeos_hybrid.h"

/*
 * Function     : NRPyEOS_hybrid_compute_rho_cold_from_P_cold()
 * Description  : Computes rho_cold; usually aliased as
 *                ghl_hybrid_compute_rho_cold_from_P_cold
 * Documentation:
*/
GHL_DEVICE
double NRPyEOS_hybrid_compute_rho_cold_from_P_cold(
      const ghl_eos_parameters *restrict eos,
      const double P_in) {

  // This code handles equations of state of the form defined
  // in Eqs 13-16 in http://arxiv.org/pdf/0802.0200.pdf
  //
  // Set up useful auxiliary variables
  const int polytropic_index = ghl_hybrid_find_polytropic_index_from_P(eos, P_in);
  const double K_ppoly       = eos->K_ppoly[polytropic_index];
  const double Gamma_ppoly   = eos->Gamma_ppoly[polytropic_index];
  double rho = pow(P_in / K_ppoly, 1.0 / Gamma_ppoly);
  NRPyEOS_hybrid_enforce_bounds__rho(eos, &rho);
  return rho;
}
