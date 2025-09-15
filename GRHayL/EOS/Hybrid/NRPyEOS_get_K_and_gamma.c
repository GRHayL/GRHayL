#include "ghl_nrpyeos_hybrid.h"

/*
 * Function     : NRPyEOS_get_K_and_Gamma()
 * Description  : For a given value of rho, find the appropriate values of
 *                Gamma_ppoly and K_ppoly by determining the appropriate index;
 *                usually aliased as ghl_hybrid_get_K_and_Gamma
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_hybrid_get_K_and_Gamma
*/

void NRPyEOS_get_K_and_Gamma(
      const ghl_eos_parameters *restrict eos,
      const double rho_in,
      double *restrict K,
      double *restrict Gamma) {

  const int polytropic_index = ghl_hybrid_find_polytropic_index(eos, rho_in);
  *K     = eos->K_ppoly[polytropic_index];
  *Gamma = eos->Gamma_ppoly[polytropic_index];
}
