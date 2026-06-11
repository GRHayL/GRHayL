#include "ghl_nrpyeos_hybrid.h"

/**
 * @ingroup hyb_eos
 * @brief Computes entropy; usually aliased as ghl_hybrid_compute_entropy_function
 *
 * @details
 * This function computes the entropy using
 *
 * \f$ S = \frac{P}{\rho^{\Gamma - 1}} \f$
 *
 * See eq. (20) in https://arxiv.org/pdf/0808.3140.pdf for more details.
 * @todo use doxygen citation
 *
 * @param[in] eos pointer to ghl_eos_parameters struct
 *
 * @param[in] rho density value
 *
 * @param[in] P pressure value
 *
 * @returns entropy
 */
double NRPyEOS_compute_entropy_function(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double P) {

  const int index    = NRPyEOS_find_polytropic_index(eos, rho);
  const double Gamma = eos->Gamma_ppoly[index];

  return P / pow(rho,Gamma-1.0);
}
