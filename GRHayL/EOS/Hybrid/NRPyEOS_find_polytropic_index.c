#include "nrpyeos_hybrid.h"

/*
 * Function     : NRPyEOS_find_polytropic_index()
 * Description  : For a given value of rho, find the appropriate polytropic index;
 *                usually aliased as ghl_hybrid_find_polytropic_index
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_hybrid_find_polytropic_index
*/

int NRPyEOS_find_polytropic_index(
      const ghl_eos_parameters *restrict eos,
      const double rho_in) {

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
  for(int j=0; j <= eos->neos-2; j++) polytropic_index += (rho_in >= eos->rho_ppoly[j]);

  return polytropic_index;
}
