#include "EOS_hybrid_header.h"

/* Function    : find_polytropic_K_and_Gamma_index()
 * Authors     : Leo Werneck & Zach Etienne
 * Description : For a given value of rho, find the
 *               appropriate values of Gamma_ppoly_tab
 *               and K_ppoly_tab by determining the appropriate
 *               index
 * Dependencies: None
 *
 * Inputs      : rho_in             - the value rho for which the polytropic
 *                                    EOS is needed
 *             : eos                - a struct containing the following
 *                                    relevant quantities:
 *             : neos               - number of polytropic EOSs used
 *             : rho_ppoly_tab      - array of rho values that determine
 *
 * Outputs     : polytropic_index   - the appropriate index for the K_ppoly_tab
 *                                    and Gamma_ppoly_tab array
 */
int find_polytropic_K_and_Gamma_index(const eos_parameters *restrict eos, double rho_in) {

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
  for(int j=0; j <= eos->neos-2; j++) polytropic_index += (rho_in >= eos->rho_ppoly_tab[j]);

  return polytropic_index;
}
