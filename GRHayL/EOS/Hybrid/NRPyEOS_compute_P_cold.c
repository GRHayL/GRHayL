#include "ghl_nrpyeos_hybrid.h"

/*
 * Function     : NRPyEOS_compute_P_cold()
 * Description  : Computes P_cold; usually aliased as
 *                ghl_hybrid_compute_P_cold_and_eps_cold
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_hybrid_compute_P_cold
*/

void NRPyEOS_compute_P_cold(
      const ghl_eos_parameters *restrict eos,
      const double rho_in,
      double *restrict P_cold_ptr) {
  // This code handles equations of state of the form defined
  // in Eqs 13-16 in http://arxiv.org/pdf/0802.0200.pdf
  //
  // Get EOS (K, Gamma)
  double K_ppoly, Gamma_ppoly;
  ghl_hybrid_get_K_and_Gamma(eos, rho_in, &K_ppoly, &Gamma_ppoly);

  // Then compute P_{cold}
  *P_cold_ptr = K_ppoly*pow(rho_in,Gamma_ppoly);
}
