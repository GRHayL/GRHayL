#include "ghl_nrpyeos_hybrid.h"

/*
 * Function     : NRPyEOS_compute_P_cold_and_eps_cold()
 * Description  : Computes P_cold and eps_cold; usually aliased as
 *                ghl_hybrid_compute_P_cold_and_eps_cold
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_hybrid_compute_P_cold_and_eps_cold
*/
GHL_DEVICE
void NRPyEOS_compute_P_cold_and_eps_cold(
      const ghl_eos_parameters *restrict eos,
      const double rho_in,
      double *restrict P_cold_ptr,
      double *restrict eps_cold_ptr) {
  // This code handles equations of state of the form defined
  // in Eqs 13-16 in http://arxiv.org/pdf/0802.0200.pdf
  //
  // Set up useful auxiliary variables
  double rho_bounded = rho_in;
  NRPyEOS_hybrid_enforce_bounds__rho(eos, &rho_bounded);

  const int polytropic_index   = ghl_hybrid_find_polytropic_index(eos, rho_bounded);
  const double K_ppoly         = eos->K_ppoly[polytropic_index];
  const double Gamma_ppoly     = eos->Gamma_ppoly[polytropic_index];
  const double eps_integ_const = eos->eps_integ_const[polytropic_index];

  // Then compute P_{cold}
  const double P_cold = K_ppoly*pow(rho_bounded, Gamma_ppoly);

  *P_cold_ptr   = P_cold;
  *eps_cold_ptr = P_cold/(rho_bounded*(Gamma_ppoly-1.0)) + eps_integ_const;
}
