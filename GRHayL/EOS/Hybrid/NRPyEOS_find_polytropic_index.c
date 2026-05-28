#include "ghl_nrpyeos_hybrid.h"

/**
 * @ingroup hyb_eos
 * @brief For a given density, find the appropriate polytropic index;
 *        usually aliased as ghl_hybrid_find_polytropic_index
 *
 * @details
 * This function finds the appropriate polytropic EOS for the input density:
 *
 * \f[
 * i = \begin{cases}
 *     0 & \text{if } \rho < \rho_0 \\
 *     1 & \text{if } \rho_0 \le \rho < \rho_1 \\
 *     2 & \text{if } \rho_1 \le \rho < \rho_2 \\
 *                 \quad \vdots & \quad \vdots \\
 *     j & \text{if } \rho_{j-1} \le \rho < \rho_j
 * \end{cases}
 * \f]
 *
 * Then, a simple way of determining the index is through the formula
 * 
 * \f$ i = (\rho \ge \rho_0) + (\rho \ge \rho_1) + \vdots + (\rho \ge \rho_{j-1})
 *
 * @param[in] eos:    pointer to ghl_eos_parameters struct
 *
 * @param[in] rho_in: density value
 *
 * @returns the index corresponding to the polytropic piece for the given density
 */
int NRPyEOS_find_polytropic_index(
      const ghl_eos_parameters *restrict eos,
      const double rho_in) {

  if(eos->neos == 1) return 0;

  int polytropic_index = 0;
  for(int j=0; j <= eos->neos-2; j++) polytropic_index += (rho_in >= eos->rho_ppoly[j]);

  return polytropic_index;
}
