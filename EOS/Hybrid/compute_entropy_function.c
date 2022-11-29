#include "EOS_hybrid_header.h"

void compute_entropy_function( const eos_parameters *restrict eos,
                               const double rho,
                               const double P,
                               double *restrict S ) {
  // This function sets the entropy funtion:
  //
  // S = P / rho^(Gamma-1)
  //
  // See eq. (20) in https://arxiv.org/pdf/0808.3140.pdf
  const int index = find_polytropic_K_and_Gamma_index(eos, rho);
  const double Gamma = eos->Gamma_ppoly_tab[index];
  // Now compute S
  *S = P / pow(rho,Gamma-1.0);
}
