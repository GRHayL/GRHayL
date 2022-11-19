/* These functions are from IGM. I need to add the author and other
   info once I figure out where they should really live. */

#include "EOS_hybrid_header.h"

int find_polytropic_K_and_Gamma_index(const eos_parameters *restrict eos, const double rho_in) {

  /* We want to find the appropriate polytropic EOS for the
   * input value rho_in. Remember that:
   *
   * if rho < rho_{0}:                P_{0} , index: 0
   * if rho >= rho_{0} but < rho_{1}: P_{1} , index: 1
   * if rho >= rho_{1} but < rho_{2}: P_{2} , index: 2
   *                      ...
   * if rho >= rho_{j-1} but < rho_{j}: P_{j} , index: j
   *
   * Then, a simple way of determining the index is through
   * the formula:
   *  ---------------------------------------------------------------------------
   * | index = (rho >= rho_{0}) + (rho >= rho_{1}) + ... + (rho >= rho_{neos-2}) |
   *  ---------------------------------------------------------------------------
   */
  if(eos->neos == 1) return 0;

  int polytropic_index = 0;
  for(int j=0; j<=eos->neos-2; j++) polytropic_index += (rho_in >= eos->rho_ppoly_tab[j]);

  return polytropic_index;
}

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
