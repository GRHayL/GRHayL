#include "nrpyeos_hybrid.h"

double NRPyEOS_compute_entropy_function(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double P) {
  // This function sets the entropy funtion:
  //
  // S = P / rho^(Gamma-1)
  //
  // See eq. (20) in https://arxiv.org/pdf/0808.3140.pdf
  const int index    = NRPyEOS_find_polytropic_index(eos, rho);
  const double Gamma = eos->Gamma_ppoly[index];
  // Now compute S
  return P / pow(rho,Gamma-1.0);
}
