#include "ghl_nrpyeos_hybrid.h"

/**
 * @ingroup hyb_eos
 * @brief For a given density, looks up the appropriate values of
 *        Gamma_ppoly and K_ppoly by determining the appropriate index;
 *        usually aliased as ghl_hybrid_get_K_and_Gamma
 *
 * @param[in] eos pointer to ghl_eos_parameters struct
 *
 * @param[in] rho_in density value
 *
 * @param[out] K polytropic constant \f$ K_j \f$ for the density's EOS segment
 *
 * @param[out] Gamma adiabatic index \f$ \Gamma_j \f$ for the density's EOS segment
 *
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
