#include "con2prim.h"

/* Function    : undensitize_conservatives()
 * Description : Computes undensitized conservatives using the metric
 *
 * Inputs      : metric         - metric_quantities struct with data for
 *                                the gridpoint of interest
 *             : cons           - conservative_quantities struct with data
 *                                for the gridpoint of interest
 *
 * Outputs     : cons_undens    - returns undensitized conservative
 *
 */

void undensitize_conservatives( const metric_quantities *restrict metric,
                                const conservative_quantities *restrict cons,
                                conservative_quantities *restrict cons_undens ) {

  // (TODO: library name)'s variables are the "densitized" versions
  // of the standard conservative variables (D,tau,S_{i}). In
  // other words, we have the relationships:
  //
  // rho_star   = sqrt(gamma) *  D
  // tilde(tau) = sqrt(gamma) * tau
  // tilde(S)_i = sqrt(gamma) * S_i
  //
  // Therefore the conversion between the two is straightfoward.

  const double psim6 = 1.0/metric->psi6;

  cons_undens->rho = cons->rho * psim6;
  cons_undens->S_x = cons->S_x * psim6;
  cons_undens->S_y = cons->S_y * psim6;
  cons_undens->S_z = cons->S_z * psim6;
  cons_undens->tau = cons->tau * psim6;
  cons_undens->Y_e = cons->Y_e * psim6;
  cons_undens->entropy = cons->entropy * psim6;
}
