#include "ghl_nrpyeos_hybrid.h"

/**
 * @ingroup hyb_eos
 * @brief For a given pressure, find the appropriate polytropic index;
 *        usually aliased as ghl_hybrid_find_polytropic_index_from_P
 *
 * @details
 * This function finds the appropriate polytropic EOS for the input pressure:
 *
 * \f[
 * i = \begin{cases}
 *     0 & \text{if } P < P_0 \\
 *     1 & \text{if } P_0 \le P < P_1 \\
 *     2 & \text{if } P_1 \le P < P_2 \\
 *                 \quad \vdots & \quad \vdots \\
 *     j & \text{if } P_{j-1} \le P < P_j
 * \end{cases}
 * \f]
 *
 * Then, a simple way of determining the index is through the formula
 * 
 * \f$ i = (P \ge P_0) + (P \ge P_1) + \vdots + (P \ge P_{j-1})
 *
 * @param[in] eos:  pointer to ghl_eos_parameters struct
 *
 * @param[in] P_in: pressure value
 *
 * @returns the index corresponding to the polytropic piece for the given pressure
 */
int NRPyEOS_find_polytropic_index_from_P(
    const ghl_eos_parameters *restrict eos,
    const double P_in) {

    if(eos->neos == 1) return 0;

    int polytropic_index = 0;
    for(int j=0; j <= eos->neos-2; j++) polytropic_index += (P_in >= eos->p_ppoly[j]);

    return polytropic_index;
}
