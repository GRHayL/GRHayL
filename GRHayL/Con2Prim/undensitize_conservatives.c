#include "ghl_con2prim.h"

/* Function     : ghl_undensitize_conservatives()
 * Description  : Computes undensitized conservatives using the metric
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_undensitize_conservatives
*/

/**
 * @ingroup Con2Prim
 * @brief Computes undensitized conservatives using the metric quantity \f$ \psi^6 \f$
 *
 * @details
 * This function computes the undensitized conservative variables from the densitized
 * conservatives using \f$ \psi^6 = \sqrt{|\gamma|} \f$ (i.e. the \f$ \psi^6 \f$ from
 * the BSSN formulation or the square root of the determinant of the ADM metric). This
 * is usually taken from ghl_metric_quantities::sqrt_detgamma, though note that this
 * instance must contain the ADM metric.
 *
 * @param[in] psi6: spacetime quantitity \f$ \psi^6 = \sqrt{|\gamma|} \f$
 *
 * @param[in] cons: pointer to ghl_conservative_quantities containing the densitized
 *                  conservative variables
 *
 * @param[in] cons_undens: pointer to ghl_conservative_quantities containing the
 *                         **undensitized** conservative variables
 *
 * @returns void
 */
void ghl_undensitize_conservatives(
      const double psi6,
      const ghl_conservative_quantities *restrict cons,
      ghl_conservative_quantities *restrict cons_undens) {

  const double psim6 = 1.0/psi6;

  cons_undens->rho     = cons->rho * psim6;
  cons_undens->SD[0]   = cons->SD[0] * psim6;
  cons_undens->SD[1]   = cons->SD[1] * psim6;
  cons_undens->SD[2]   = cons->SD[2] * psim6;
  cons_undens->tau     = cons->tau * psim6;
  cons_undens->Y_e     = cons->Y_e * psim6;
  cons_undens->entropy = cons->entropy * psim6;
}
