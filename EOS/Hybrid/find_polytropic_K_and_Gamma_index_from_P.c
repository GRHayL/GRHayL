#include "EOS_hybrid_header.h"

/* Function    : find_polytropic_K_and_Gamma_index_from_P()
 * Authors     : Leo Werneck
 * Description : For a given value of , find the
 *               appropriate values of Gamma_ppoly_tab
 *               and K_ppoly_tab by determining the appropriate
 *               index
 * Dependencies: initialize_igm_eos_parameters_from_input()
 *             : cctk_parameters.h (FIXME)
 *
 * Inputs      : eos               - a struct containing the following
 *                                   relevant quantities:
 *                 : neos          - number of polytropic EOSs used
 *                 : rho_ppoly_tab - array of rho values that determine
 *             : P_in              - Input pressure
 *
 * Outputs     : index             - the appropriate index for the K_ppoly_tab
 *                                   and Gamma_ppoly_tab array
 */
int find_polytropic_K_and_Gamma_index_from_P(const igm_eos_parameters eos, const CCTK_REAL P_in) {

  if(eos.neos == 1) return 0;

  int polytropic_index = 0;
  for(int j=0; j<=eos.neos-2; j++) {
    // This function is just slightly more involved than the previous function.
    // Instead of comparing rho_in against the values of rho that separate
    // one polytropic index from the next, we compute P(rho) at the separation
    // values and compare P_in against them. This is guaranteed to work because
    // P(rho) increases monotonically with rho.
    const CCTK_REAL P_local = eos.K_ppoly_tab[j] * pow( eos.rho_ppoly_tab[j], eos.Gamma_ppoly_tab[j] );
  //Original version: for(int j=0; j<=eos->neos-2; j++) polytropic_index += (rho_in >= eos->rho_ppoly_tab[j]);
    polytropic_index += (P_in >= P_local);
  }

  return polytropic_index;
}
