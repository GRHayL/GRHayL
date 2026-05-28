#include "ghl_nrpyeos_hybrid.h"

/**
 * @ingroup hyb_eos
 * @brief Computes cold pressure; usually aliased as ghl_hybrid_compute_P_cold
 *
 * @details
 * This code handles equations of state of the form defined
 * in Eqs 13-16 in http://arxiv.org/pdf/0802.0200.pdf
 * @todo use doxygen citation
 *
 * After finding the polytropic piece for the given density, the
 * function evaluates
 *
 * \f$ P_\mathrm{cold} = K_i \rho^{\Gamma_i} \f$
 *
 * @param[in] eos:         pointer to ghl_eos_parameters struct
 *
 * @param[in] rho_in:      density value
 *
 * @param[out] P_cold_ptr: pointer to ghl_eos_parameters struct
 *
 * @returns void
 */
void NRPyEOS_compute_P_cold(
      const ghl_eos_parameters *restrict eos,
      const double rho_in,
      double *restrict P_cold_ptr) {

  // Get EOS (K, Gamma)
  double K_ppoly, Gamma_ppoly;
  ghl_hybrid_get_K_and_Gamma(eos, rho_in, &K_ppoly, &Gamma_ppoly);

  // Then compute P_{cold}
  *P_cold_ptr = K_ppoly*pow(rho_in,Gamma_ppoly);
}
