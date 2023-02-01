#include "NRPyEOS_Hybrid.h"

/* Function    : get_K_and_Gamma()
 * Description : For a given value of rho, find the
 *               appropriate values of Gamma_ppoly
 *               and K_ppoly by determining the
 *               appropriate index
 * Dependencies: None
 *
 * Inputs      : rho_in         - the value rho for which the polytropic
 *                                EOS is needed
 *             : eos            - an initialized eos_parameters struct
 *                                with data for the EOS of the simulation
 *
 * Outputs     : K              - the value of polytropic K for rho_in
 * Outputs     : Gamma          - the value of polytropic gamma for rho_in
 */
void NRPyEOS_get_K_and_Gamma(
      const eos_parameters *restrict eos,
      const double rho_in,
      double *restrict K,
      double *restrict Gamma) {

  int polytropic_index = NRPyEOS_find_polytropic_index(eos, rho_in);
  *K     = eos->K_ppoly[polytropic_index];
  *Gamma = eos->Gamma_ppoly[polytropic_index];
}
