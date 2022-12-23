#include "GRHayL.h"

void GRHayL_compute_P_cold(
           const eos_parameters *restrict eos,
           const double rho_in,
           double *restrict P_cold_ptr) {
  // This code handles equations of state of the form defined
  // in Eqs 13-16 in http://arxiv.org/pdf/0802.0200.pdf
  //
  // Get EOS (K, Gamma)
  double K_ppoly, Gamma_ppoly;
  (*eos->hybrid_get_K_and_Gamma)(eos, rho_in, &K_ppoly, &Gamma_ppoly);

  // Then compute P_{cold}
  *P_cold_ptr = K_ppoly*pow(rho_in,Gamma_ppoly);
}
